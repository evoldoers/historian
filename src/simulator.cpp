#include "simulator.h"
#include "util.h"
#include "logger.h"

AlignPath Simulator::simulateGapsByGillespie (random_engine& generator, const RateModel& model, SeqIdx parentLength, double time, AlignRowIndex parentRowIndex, AlignRowIndex childRowIndex) {
  LogThisAt(4,"Simulating indels in " << parentLength << "-residue sequence for branch length " << time << endl);
  vguard<int> child2parent (parentLength);
  iota (child2parent.begin(), child2parent.end(), 0);
  double t = 0;
  uniform_real_distribution<double> uniform (0, 1);
  geometric_distribution<int> insLenDist (1. - model.insExtProb);
  geometric_distribution<int> delLenDist (1. - model.delExtProb);
  while (true) {
    const double totalInsRate = model.insRate * (child2parent.size() + 1);
    const double totalDelRate = model.delRate * child2parent.size();
    const double totalIndelRate = totalInsRate + totalDelRate;
    const double wait = -log(uniform(generator)) / totalIndelRate;
    LogThisAt(5,"Current time is " << t << ", wait for next event is " << wait << endl);
    t += wait;
    if (t > time) {
      LogThisAt(5,"Passed time limit; stopping" << endl);
      break;
    }
    const double r = uniform(generator) * totalIndelRate;
    const bool isInsertion = r < totalInsRate;
    if (isInsertion) {
      const SeqIdx insPos = (SeqIdx) (r / model.insRate);
      const SeqIdx insLen = 1 + insLenDist (generator);
      LogThisAt(5,"Insertion of " << plural(insLen,"residue") << " before position #" << insPos << endl);
      child2parent.insert (child2parent.begin() + insPos, insLen, -1);
    } else {
      const double rDel = r - totalInsRate;
      const SeqIdx delPos = (SeqIdx) (rDel / model.delRate);
      const SeqIdx delLen = 1 + delLenDist (generator);
      const SeqIdx delEnd = min (delPos + delLen, (SeqIdx) child2parent.size());
      LogThisAt(5,"Deletion of " << plural(delEnd-delPos,"residue") << " starting at position #" << delPos << endl);
      child2parent.erase (child2parent.begin() + delPos, child2parent.begin() + delEnd);
    }
  }
  LogThisAt(6,"Child-to-parent map: " << to_string_join(child2parent) << endl);
  AlignRowPath parentPath, childPath;
  SeqIdx parentPos = 0;
  for (SeqIdx childPos = 0; childPos < child2parent.size(); ++childPos) {
    if (child2parent[childPos] < 0) {
      parentPath.push_back (false);
      childPath.push_back (true);
    } else {
      while (parentPos < (SeqIdx) child2parent[childPos]) {
	parentPath.push_back (true);
	childPath.push_back (false);
	++parentPos;
      }
      parentPath.push_back (true);
      childPath.push_back (true);
      ++parentPos;
    }
  }
  while (parentPos < parentLength) {
    parentPath.push_back (true);
    childPath.push_back (false);
    ++parentPos;
  }
  AlignPath path;
  path[parentRowIndex] = parentPath;
  path[childRowIndex] = childPath;
  LogThisAt(6,"Child-parent alignment:\n" << alignPathString(path) << endl);
  return path;
}

vguard<FastSeq> Simulator::simulateSubsByMatExp (random_engine& generator, const RateModel& model, const Tree& tree, const AlignPath& path) {
  const AlignColIndex cols = alignPathColumns(path);
  const AlignRowIndex rows = tree.nodes();
  vguard<FastSeq> gapped (rows);
  vguard<vguard<AlphTok> > tok (rows, cols);
  vguard<vguard<int> > component (rows, cols);
  discrete_distribution<int> cptDist (model.cptWeight.begin(), model.cptWeight.end());
  vguard<discrete_distribution<AlphTok> > cptInsDist (model.components());
  vguard<vguard<vguard<discrete_distribution<AlphTok> > > > nodeCptCondSubDist (tree.nodes(),
										vguard<vguard<discrete_distribution<AlphTok> > > (model.components(),
																  vguard<discrete_distribution<AlphTok> > (model.alphabetSize())));
  for (size_t c = 0; c < model.components(); ++c) {
    const auto insProb = gsl_vector_to_stl (model.insProb[c]);
    cptInsDist[c] = discrete_distribution<AlphTok> (insProb.begin(), insProb.end());
  }
  for (auto node: tree.preorderSort()) {
    gapped[node].name = tree.seqName(node);
    gapped[node].seq = string (cols, Alignment::gapChar);
    gapped[node].qual = string (cols, Alignment::gapChar);
    vguard<gsl_matrix*> subMat = model.getSubProbMatrix (tree.branchLength (node));
    for (size_t c = 0; c < model.components(); ++c) {
      vguard<vguard<double> > cptSubMat = gsl_matrix_to_stl (subMat[c]);
      for (AlphTok i = 0; i < model.alphabetSize(); ++i)
	nodeCptCondSubDist[node][c][i] = discrete_distribution<AlphTok> (cptSubMat[i].begin(), cptSubMat[i].end());
      gsl_matrix_free (subMat[c]);
    }
  }
  for (auto node: tree.preorderSort())
    for (AlignRowIndex col = 0; col < cols; ++col) {
      if (path.at(node)[col]) {
	auto parent = tree.parentNode (node);
	const bool isInsertion = parent < 0 || !path.at(parent)[col];
	int cpt;
	if (isInsertion) {
	  cpt = cptDist (generator);
	  tok[node][col] = cptInsDist[cpt] (generator);
	} else {
	  cpt = component[parent][col];
	  tok[node][col] = nodeCptCondSubDist[node][cpt][tok[parent][col]] (generator);
	}
	component[node][col] = cpt;
	gapped[node].seq[col] = model.alphabet[tok[node][col]];
	gapped[node].qual[col] = cpt < 10 ? ('0' + cpt) : ('A' + cpt - 10);
      }
    }
  return gapped;
}

Stockholm Simulator::simulateTree (random_engine& generator, const RateModel& model, const Tree& tree, SeqIdx rootLength) {
  vguard<AlignPath> branchPaths;
  vguard<SeqIdx> nodeLen (tree.nodes());
  vguard<FastSeq> wild (tree.nodes());
  nodeLen[tree.root()] = rootLength;
  for (TreeNodeIndex node = tree.root() - 1; node >= 0; --node) {
    const auto parent = tree.parentNode(node);
    LogThisAt(4,"Simulating branch from " << tree.seqName(parent) << " to " << tree.seqName(node) << endl);
    const auto parentLen = nodeLen[parent];
    const auto branchPath = simulateGapsByGillespie (generator, model, parentLen, tree.branchLength(node), parent, node);
    nodeLen[node] = alignPathResiduesInRow (branchPath.at(node));
    Assert (parentLen == alignPathResiduesInRow (branchPath.at(parent)), "Parent length mismatch");
    branchPaths.push_back (branchPath);
  }
  for (TreeNodeIndex node = 0; node < tree.nodes(); ++node) {
    wild[node].name = tree.seqName(node);
    wild[node].seq = string (rootLength, Alignment::wildcardChar);
  }
  const AlignPath path = alignPathMerge (branchPaths);
  const vguard<FastSeq> gapped = simulateSubsByMatExp (generator, model, tree, path);
  Stockholm stock (gapped, tree);
  stock.assertFlush();
  if (model.components() > 1)
    for (TreeNodeIndex node = 0; node < tree.nodes(); ++node)
      stock.gr[string(SimulatorComponentTag)][stock.gapped[node].name] = stock.gapped[node].qual;
  for (auto& fs: stock.gapped)
    fs.qual.clear();
  return stock;
}
