#include <fstream>
#include <random>
#include "recon.h"
#include "util.h"
#include "forward.h"
#include "jsonutil.h"
#include "logger.h"
#include "span.h"
#include "amino.h"

Reconstructor::Reconstructor()
  : profileSamples (DefaultProfileSamples),
    profileNodeLimit (0),
    rndSeed (ForwardMatrix::random_engine::default_seed),
    maxDistanceFromGuide (DefaultMaxDistanceFromGuide),
    includeBestTraceInProfile (true),
    keepGapsOpen (false),
    usePosteriorsForProfile (true),
    reconstructRoot (true),
    accumulateCounts (false),
    predictAncestralSequence (false),
    minPostProb (DefaultProfilePostProb)
{ }

bool Reconstructor::parseReconArgs (deque<string>& argvec) {
  if (argvec.size()) {
    const string& arg = argvec[0];
    if (arg == "-ancseq") {
      predictAncestralSequence = true;
      argvec.pop_front();
      return true;
    }
  }
  return parsePostArgs (argvec);
}

bool Reconstructor::parsePostArgs (deque<string>& argvec) {
  if (argvec.size()) {
    const string& arg = argvec[0];
    if (arg == "-seqs") {
      Require (argvec.size() > 1, "%s must have an argument", arg.c_str());
      Require (guideFilename.size() == 0, "Can't specify both -guide and -seqs");
      seqsFilename = argvec[1];
      argvec.pop_front();
      argvec.pop_front();
      return true;

    } else if (arg == "-guide") {
      Require (argvec.size() > 1, "%s must have an argument", arg.c_str());
      Require (seqsFilename.size() == 0, "Can't specify both -seqs and -guide");
      guideFilename = argvec[1];
      argvec.pop_front();
      argvec.pop_front();
      return true;

    } else if (arg == "-savetree") {
      Require (argvec.size() > 1, "%s must have an argument", arg.c_str());
      treeSaveFilename = argvec[1];
      argvec.pop_front();
      argvec.pop_front();
      return true;

    } else if (arg == "-saveseqs") {
      Require (argvec.size() > 1, "%s must have an argument", arg.c_str());
      seqsSaveFilename = argvec[1];
      argvec.pop_front();
      argvec.pop_front();
      return true;

    } else if (arg == "-saveguide") {
      Require (argvec.size() > 1, "%s must have an argument", arg.c_str());
      guideSaveFilename = argvec[1];
      argvec.pop_front();
      argvec.pop_front();
      return true;

    } else if (arg == "-band") {
      Require (argvec.size() > 1, "%s must have an argument", arg.c_str());
      maxDistanceFromGuide = atoi (argvec[1].c_str());
      argvec.pop_front();
      argvec.pop_front();
      return true;

    } else if (arg == "-noband") {
      maxDistanceFromGuide = -1;
      argvec.pop_front();
      return true;

    } else if (arg == "-samples") {
      Require (argvec.size() > 1, "%s must have an argument", arg.c_str());
      profileSamples = atoi (argvec[1].c_str());
      argvec.pop_front();
      argvec.pop_front();
      usePosteriorsForProfile = false;
      return true;

    } else if (arg == "-minpost") {
      Require (argvec.size() > 1, "%s must have an argument", arg.c_str());
      minPostProb = atof (argvec[1].c_str());
      argvec.pop_front();
      argvec.pop_front();
      usePosteriorsForProfile = true;
      return true;
      
    } else if (arg == "-states") {
      Require (argvec.size() > 1, "%s must have an argument", arg.c_str());
      profileNodeLimit = atoi (argvec[1].c_str());
      argvec.pop_front();
      argvec.pop_front();
      return true;

    } else if (arg == "-nobest") {
      includeBestTraceInProfile = false;
      argvec.pop_front();
      return true;

    } else if (arg == "-keepgapsopen") {
      keepGapsOpen = true;
      argvec.pop_front();
      return true;

    } else if (arg == "-seed") {
      Require (argvec.size() > 1, "%s must have an argument", arg.c_str());
      rndSeed = atoi (argvec[1].c_str());
      argvec.pop_front();
      argvec.pop_front();
      return true;
    }
  }

  return diagEnvParams.parseDiagEnvParams (argvec)
    || parseTreeArgs (argvec)
    || parseModelArgs (argvec);
}

bool Reconstructor::parseCountArgs (deque<string>& argvec) {
  if (argvec.size()) {
    const string& arg = argvec[0];
    if (arg == "-recon") {
      Require (argvec.size() > 1, "%s must have an argument", arg.c_str());
      Require (reconFilename.size() == 0, "Can't specify both -guide and -seqs");
      reconFilename = argvec[1];
      argvec.pop_front();
      argvec.pop_front();
      return true;
    }
  }

  return parseTreeArgs (argvec)
    || parseModelArgs (argvec);
}

bool Reconstructor::parseTreeArgs (deque<string>& argvec) {
  if (argvec.size()) {
    const string& arg = argvec[0];
    if (arg == "-tree") {
      Require (argvec.size() > 1, "%s must have an argument", arg.c_str());
      treeFilename = argvec[1];
      argvec.pop_front();
      argvec.pop_front();
      return true;
    }
  }

  return false;
}

bool Reconstructor::parseModelArgs (deque<string>& argvec) {
  if (argvec.size()) {
    const string& arg = argvec[0];
    if (arg == "-model") {
      Require (argvec.size() > 1, "%s must have an argument", arg.c_str());
      modelFilename = argvec[1];
      argvec.pop_front();
      argvec.pop_front();
      return true;

    } else if (arg == "-savemodel") {
      Require (argvec.size() > 1, "%s must have an argument", arg.c_str());
      modelSaveFilename = argvec[1];
      argvec.pop_front();
      argvec.pop_front();
      return true;
    }
  }

  return false;
}

bool Reconstructor::parseSumArgs (deque<string>& argvec) {
  if (argvec.size()) {
    const string& arg = argvec[0];
    if (arg == "-counts") {
      Require (argvec.size() > 1, "%s must have an argument", arg.c_str());
      countFilenames.push_back (argvec[1]);
      argvec.pop_front();
      argvec.pop_front();
      return true;
    }
  }

  return false;
}

void Reconstructor::loadModel() {
  if (modelFilename.size()) {
    LogThisAt(1,"Loading model from " << modelFilename << endl);
    ifstream modelFile (modelFilename);
    ParsedJson pj (modelFile);
    model.read (pj.value);
  } else {
    LogThisAt(1,"Using default amino acid model" << endl);
    model = defaultAminoModel();
  }
}

void Reconstructor::loadTree() {
  Require (treeFilename.size() > 0, "Must specify a tree");
  LogThisAt(1,"Loading tree from " << treeFilename << endl);
  ifstream treeFile (treeFilename);
  dataset.tree.parse (JsonUtil::readStringFromStream (treeFile));
}

void Reconstructor::buildTree() {
  LogThisAt(1,"Building neighbor-joining tree" << endl);
  auto dist = model.distanceMatrix (dataset.gappedGuide);
  dataset.tree.buildByNeighborJoining (dataset.gappedGuide, dist);
}

void Reconstructor::seedGenerator() {
  generator = ForwardMatrix::newRNG();
  generator.seed (rndSeed);
}

void Reconstructor::loadSeqs() {
  Require (seqsFilename.size() > 0 || guideFilename.size() > 0, "Must specify sequences");

  loadModel();
  seedGenerator();
  
  if (seqsFilename.size()) {
    LogThisAt(1,"Loading sequences from " << seqsFilename << endl);
    dataset.seqs = readFastSeqs (seqsFilename.c_str());
    if (maxDistanceFromGuide < 0 && treeFilename.size())
      LogThisAt(1,"Don't need guide alignment: banding is turned off and tree is supplied" << endl);
    else {
      LogThisAt(1,"Building guide alignment" << endl);
      AlignGraph ag (dataset.seqs, model, 1, diagEnvParams, generator);
      Alignment align = ag.mstAlign();
      dataset.guide = align.path;
      dataset.gappedGuide = align.gapped();
    }

  } else {
    LogThisAt(1,"Loading guide alignment from " << guideFilename << endl);
    dataset.gappedGuide = readFastSeqs (guideFilename.c_str());
    const Alignment align (dataset.gappedGuide);
    dataset.guide = align.path;
    dataset.seqs = align.ungapped;
  }

  if (treeFilename.size())
    loadTree();
  else
    buildTree();

  dataset.buildReconIndices();

  if (modelSaveFilename.size()) {
    ofstream modelFile (modelSaveFilename);
    model.write (modelFile);
  }

  if (seqsSaveFilename.size()) {
    ofstream seqsFile (seqsSaveFilename);
    writeFastaSeqs (seqsFile, dataset.seqs);
  }

  if (guideSaveFilename.size()) {
    if (dataset.gappedGuide.empty())
      Warn("No guide alignment to save");
    else {
      ofstream guideFile (guideSaveFilename);
      writeFastaSeqs (guideFile, dataset.gappedGuide);
    }
  }

  if (treeSaveFilename.size()) {
    ofstream treeFile (treeSaveFilename);
    treeFile << dataset.tree.toString() << endl;
  }
}

void Reconstructor::Dataset::buildReconIndices() {
  for (size_t n = 0; n < seqs.size(); ++n) {
    Assert (seqIndex.find (seqs[n].name) == seqIndex.end(), "Duplicate sequence name %s", seqs[n].name.c_str());
    seqIndex[seqs[n].name] = n;
  }

  for (TreeNodeIndex node = 0; node < tree.nodes(); ++node)
    if (!tree.isLeaf(node))
      Assert (tree.nChildren(node) == 2, "Tree is not binary: node %d has %s\nSubtree rooted at %d: %s\n", node, plural(tree.nChildren(node),"child","children").c_str(), node, tree.toString(node).c_str());

  AlignPath reorderedGuide;
  for (TreeNodeIndex node = 0; node < tree.nodes(); ++node)
    if (tree.isLeaf(node)) {
      Assert (tree.nodeName(node).length() > 0, "Leaf node %d is unnamed", node);
      Assert (seqIndex.find (tree.nodeName(node)) != seqIndex.end(), "Can't find sequence for leaf node %s", tree.nodeName(node).c_str());
      const size_t seqidx = seqIndex[tree.nodeName(node)];
      nodeToSeqIndex[node] = seqidx;

      if (!guide.empty())
	reorderedGuide[node] = guide.at(seqidx);

      closestLeaf.push_back (node);
      closestLeafDistance.push_back (0);

      rowName.push_back (tree.seqName(node));

    } else {
      int cl = -1;
      double dcl = 0;
      for (size_t nc = 0; nc < tree.nChildren(node); ++nc) {
	const TreeNodeIndex c = tree.getChild(node,nc);
	const double dc = closestLeafDistance[c] + tree.branchLength(c);
	if (nc == 0 || dc < dcl) {
	  cl = closestLeaf[c];
	  dcl = dc;
	}
      }
      closestLeaf.push_back (cl);
      closestLeafDistance.push_back (dcl);

      rowName.push_back (tree.seqName(node));
    }

  swap (guide, reorderedGuide);
}

void Reconstructor::reconstruct() {
  LogThisAt(1,"Starting reconstruction on " << dataset.tree.nodes() << "-node tree" << endl);

  seedGenerator();  // re-seed generator, in case it was used during prealignment

  gsl_vector* rootProb = model.insProb;
  LogProb lpFinalFwd = -numeric_limits<double>::infinity(), lpFinalTrace = -numeric_limits<double>::infinity();
  const ForwardMatrix::ProfilingStrategy strategy =
    (ForwardMatrix::ProfilingStrategy) (ForwardMatrix::CollapseChains
					| (keepGapsOpen ? ForwardMatrix::KeepGapsOpen : ForwardMatrix::DontKeepGapsOpen)
					| (accumulateCounts ? ForwardMatrix::CountEvents : ForwardMatrix::DontCountEvents)
					| (includeBestTraceInProfile ? ForwardMatrix::IncludeBestTrace : ForwardMatrix::DontIncludeBestTrace));

  SumProduct* sumProd = NULL;
  if (accumulateCounts)
    sumProd = new SumProduct (model, dataset.tree);

  AlignPath path;
  map<int,Profile> prof;
  for (TreeNodeIndex node = 0; node < dataset.tree.nodes(); ++node) {
    if (dataset.tree.isLeaf(node))
      prof[node] = Profile (model.alphabet, dataset.seqs[dataset.nodeToSeqIndex[node]], node);
    else {
      const int lChildNode = dataset.tree.getChild(node,0);
      const int rChildNode = dataset.tree.getChild(node,1);
      const Profile& lProf = prof[lChildNode];
      const Profile& rProf = prof[rChildNode];
      ProbModel lProbs (model, dataset.tree.branchLength(lChildNode));
      ProbModel rProbs (model, dataset.tree.branchLength(rChildNode));
      PairHMM hmm (lProbs, rProbs, rootProb);

      LogThisAt(2,"Aligning " << lProf.name << " (" << plural(lProf.state.size(),"state") << ", " << plural(lProf.trans.size(),"transition") << ") and " << rProf.name << " (" << plural(rProf.state.size(),"state") << ", " << plural(rProf.trans.size(),"transition") << ")" << endl);

      ForwardMatrix forward (lProf, rProf, hmm, node, dataset.guide.empty() ? GuideAlignmentEnvelope() : GuideAlignmentEnvelope (dataset.guide, dataset.closestLeaf[lChildNode], dataset.closestLeaf[rChildNode], maxDistanceFromGuide), sumProd);

      BackwardMatrix *backward = NULL;
      if ((accumulateCounts && node == dataset.tree.root()) || (usePosteriorsForProfile && node != dataset.tree.root()))
	backward = new BackwardMatrix (forward, minPostProb);

      Profile& nodeProf = prof[node];
      if (node == dataset.tree.root()) {
	if (reconstructRoot) {
	  path = forward.bestAlignPath();
	  nodeProf = forward.bestProfile();
	}
      } else if (usePosteriorsForProfile)
	nodeProf = backward->buildProfile (profileNodeLimit, strategy);
      else
	nodeProf = forward.sampleProfile (generator, profileSamples, profileNodeLimit, strategy);

      if (accumulateCounts && node == dataset.tree.root())
	dataset.eigenCounts = backward->getCounts();
      
      if (backward)
	delete backward;
      
      if (node == dataset.tree.root())
	lpFinalFwd = forward.lpEnd;

      if (nodeProf.size()) {
        const LogProb lpTrace = nodeProf.calcSumPathAbsorbProbs (log_gsl_vector(rootProb), NULL);
        LogThisAt(3,"Forward log-likelihood is " << forward.lpEnd << ", profile log-likelihood is " << lpTrace << " with " << nodeProf.size() << " states" << endl);

	if (node == dataset.tree.root())
	  lpFinalTrace = lpTrace;

	LogThisAt(5,nodeProf.toJson());
      }
    }
  }

  LogThisAt(1,"Final Forward log-likelihood is " << lpFinalFwd << (reconstructRoot ? (string(", final alignment log-likelihood is ") + to_string(lpFinalTrace)) : string()) << endl);

  if (reconstructRoot) {
    vguard<FastSeq> ungapped (dataset.tree.nodes());
    for (TreeNodeIndex node = 0; node < dataset.tree.nodes(); ++node) {
      if (dataset.tree.isLeaf(node)) 
	ungapped[node] = dataset.seqs[dataset.seqIndex.at(dataset.rowName[node])];
      else {
	ungapped[node].seq = string (alignPathResiduesInRow(path.at(node)), '*');
	ungapped[node].name = dataset.rowName[node];
      }
    }

    dataset.reconstruction = Alignment (ungapped, path);
    dataset.gappedRecon = dataset.reconstruction.gapped();

    if (predictAncestralSequence) {
      AlignColSumProduct colSumProd (model, dataset.tree, dataset.gappedRecon);
      while (!colSumProd.alignmentDone()) {
	colSumProd.fillUp();
	colSumProd.fillDown();
	colSumProd.appendAncestralReconstructedColumn (dataset.ancestral);
	colSumProd.nextColumn();
      }
    }
  }

  if (accumulateCounts)
    dataset.eventCounts = dataset.eigenCounts.transform (model);
  
  if (sumProd)
    delete sumProd;
}

void Reconstructor::writeRecon (ostream& out) const {
  writeFastaSeqs (cout, predictAncestralSequence ? dataset.ancestral : dataset.gappedRecon);
}

void Reconstructor::writeCounts (ostream& out) const {
  dataset.eventCounts.writeJson (out);
}

void Reconstructor::writeModel (ostream& out) const {
  model.write (out);
}

void Reconstructor::loadRecon() {
  loadModel();
  loadTree();

  Require (reconFilename.size() > 0, "Must specify a reconstruction file");
  LogThisAt(1,"Loading reconstruction from " << reconFilename << endl);
  dataset.gappedRecon = readFastSeqs (reconFilename.c_str());
  dataset.tree.reorder (dataset.gappedRecon);

  dataset.reconstruction = Alignment (dataset.gappedRecon);
}

void Reconstructor::loadCounts() {
  for (size_t n = 0; n < countFilenames.size(); ++n) {
    ifstream in (countFilenames[n]);
    ParsedJson pj (in);
    EventCounts c;
    c.read (pj.value);
    if (n == 0)
      dataset.eventCounts = c;
    else
      dataset.eventCounts += c;
  }
}

void Reconstructor::count() {
  dataset.eigenCounts = EigenCounts (model.alphabetSize());
  dataset.eigenCounts.accumulateCounts (model, dataset.reconstruction, dataset.tree);
  dataset.eventCounts = dataset.eigenCounts.transform (model);
}

void Reconstructor::fit() {
  dataset.eventCounts.optimize (model);
}
