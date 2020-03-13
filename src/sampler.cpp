#include <gsl/gsl_math.h>
#include "sampler.h"
#include "recon.h"
#include "util.h"

#define SAMPLER_EPSILON 1e-3
#define SAMPLER_NEAR_EQ(X,Y) (gsl_fcmp (X, Y, SAMPLER_EPSILON) == 0)

double SimpleTreePrior::coalescenceRate (int lineages) const {
  return (((double) lineages * (lineages-1)) / 2) / (double) populationSize;
}

LogProb SimpleTreePrior::treeLogLikelihood (const Tree& tree) const {
  tree.assertBinary();
  const auto distanceFromRoot = tree.distanceFromRoot();
  const auto nodesByDistanceFromRoot = orderedIndices (distanceFromRoot);
  size_t lineages = 0;
  LogProb lp = 0;
  double lastEventTime;
  for (auto iter = nodesByDistanceFromRoot.rbegin(); iter != nodesByDistanceFromRoot.rend(); ++iter) {
    const double eventTime = distanceFromRoot[*iter];
    if (lineages > 1)
      lp -= coalescenceRate(lineages) * (lastEventTime - eventTime);
    lastEventTime = eventTime;
    if (tree.isLeaf(*iter))
      ++lineages;
    else
      --lineages;
  }
  return lp;
}

Sampler::History Sampler::History::reorder (const vguard<TreeNodeIndex>& newOrder) const {
  LogThisAt(6,"Reordering nodes to maintain preorder sort (" << to_string_join(newOrder) << ")" << endl);
  History newHistory;
  newHistory.tree = tree.reorderNodes (newOrder);
  newHistory.gapped.reserve (gapped.size());
  for (auto n : newOrder)
    newHistory.gapped.push_back (gapped[n]);
  return newHistory;
}

void Sampler::History::assertNamesMatch() const {
  tree.assertAllNodesNamed();
  tree.assertNodesMatchSeqs (gapped);
}

TreeNodeIndex Sampler::randomInternalNode (const Tree& tree, random_engine& generator) {
  vguard<TreeNodeIndex> intNodes;
  intNodes.reserve (tree.nodes() / 2);
  for (TreeNodeIndex n = 0; n < tree.nodes(); ++n)
    if (!tree.isLeaf(n))
      intNodes.push_back (n);
  return random_element (intNodes, generator);
}

TreeNodeIndex Sampler::randomChildNode (const Tree& tree, random_engine& generator) {
  Assert (tree.hasChildren(), "No child nodes in tree");
  std::uniform_int_distribution<TreeNodeIndex> distribution (0, tree.nodes() - 2);
  return distribution (generator);
}

TreeNodeIndex Sampler::randomGrandchildNode (const Tree& tree, random_engine& generator) {
  Assert (tree.hasGrandchildren(), "No grandchild nodes found");
  vguard<TreeNodeIndex> grandkids;
  for (TreeNodeIndex n = 0; n < tree.root(); ++n)
    if (tree.parentNode(n) != tree.root())
      grandkids.push_back (n);
  return random_element (grandkids, generator);
}

vguard<TreeNodeIndex> Sampler::contemporaneousNodes (const Tree& tree, const vguard<TreeBranchLength>& dist, TreeNodeIndex node) {
  vguard<TreeNodeIndex> contemps;
  const TreeNodeIndex parent = tree.parentNode(node);
  Assert (parent >= 0, "Parent node not found");
  Assert (tree.parentNode(parent) >= 0, "Grandparent node not found");
  const TreeBranchLength distParent = dist[parent];
  for (TreeNodeIndex n = 0; n < tree.root(); ++n) {
    const TreeNodeIndex p = tree.parentNode(n);
    if (p != parent && dist[p] < distParent && dist[n] > distParent)
      contemps.push_back (n);
  }
  const auto ndist = tree.distanceFrom (node);
  sortIndices (contemps, ndist);
  return contemps;
}

vguard<double> Sampler::nodeListWeights (size_t n) {
  vguard<double> w (n);
  double norm = 0, wi = 1, mul = 1.5;
  for (size_t i = 0; i < n; ++i) {
    w[i] = wi;
    norm += wi;
    wi /= mul;
  }
  for (auto& weight: w)
    weight /= norm;
  return w;
}

AlignRowIndex Sampler::guideRow (const Tree& tree, TreeNodeIndex node) const {
  return guideRowByName.at (tree.nodeName (node));
}

GuideAlignmentEnvelope Sampler::makeGuide (const Tree& tree, TreeNodeIndex leaf1, TreeNodeIndex leaf2, const AlignPath& path, TreeNodeIndex node1, TreeNodeIndex node2) const {
  return GuideAlignmentEnvelope (useFixedGuide ? guide.path : path, useFixedGuide ? guideRow(tree,leaf1) : node1, useFixedGuide ? guideRow(tree,leaf2) : node2, maxDistanceFromGuide);
}

vguard<SeqIdx> Sampler::guideSeqPos (const AlignPath& path, AlignRowIndex row, AlignRowIndex fixedGuideRow) const {
  return guideSeqPos (path, row, row, fixedGuideRow);
}

vguard<SeqIdx> Sampler::guideSeqPos (const AlignPath& path, AlignRowIndex row, AlignRowIndex variableGuideRow, AlignRowIndex fixedGuideRow) const {
  const AlignRowIndex guideRow = useFixedGuide ? fixedGuideRow : variableGuideRow;
  return TreeAlignFuncs::getGuideSeqPos (path, row, guideRow);
}

vguard<SeqIdx> TreeAlignFuncs::getGuideSeqPos (const AlignPath& path, AlignRowIndex row, AlignRowIndex guideRow) {
  vguard<SeqIdx> guidePos;
  LogThisAt(9,"Mapping sequence coordinates between node #" << row << " and guide row #" << guideRow << endl);
  const auto cols = alignPathColumns (path);
  guidePos.reserve (cols);
  const AlignRowPath& rowPath = path.at(row);
  const AlignRowPath& guideRowPath = path.at(guideRow);
  SeqIdx pos = 0;
  guidePos.push_back (0);
  for (AlignColIndex col = 0; col < cols; ++col) {
    if (guideRowPath[col])
      ++pos;
    if (rowPath[col])
      guidePos.push_back (pos);
  }
  return guidePos;
}

AlignPath TreeAlignFuncs::cladePath (const AlignPath& path, const Tree& tree, TreeNodeIndex cladeRoot, TreeNodeIndex cladeRootParent, TreeNodeIndex exclude) {
  AlignPath p;
  const vguard<TreeNodeIndex> rerootedParent = tree.rerootedParent (cladeRootParent);
  vguard<bool> childrenIncluded (tree.nodes(), false);
  childrenIncluded[cladeRootParent] = true;
  for (auto n : tree.rerootedPreorderSort(cladeRoot,cladeRootParent)) {
    if (n != exclude && childrenIncluded[rerootedParent[n]]) {
      p[n] = path.at(n);
      childrenIncluded[n] = true;
    }
  }
  return alignPathRemoveEmptyColumns (p);
}

AlignPath TreeAlignFuncs::pairPath (const AlignPath& path, TreeNodeIndex node1, TreeNodeIndex node2) {
  AlignPath p;
  size_t nDel = 0;
  const AlignColIndex cols = alignPathColumns (path);
  const AlignRowPath& row1 = path.at (node1);
  const AlignRowPath& row2 = path.at (node2);
  AlignRowPath& r1 = p[node1];
  AlignRowPath& r2 = p[node2];
  for (AlignColIndex col = 0; col < cols; ++col) {
    const bool c1 = row1[col];
    const bool c2 = row2[col];
    if (!(c1 || c2))
      continue;
    const ProbModel::State state = ProbModel::getState (c1, c2);
    switch (state) {
    case ProbModel::Match:
      while (nDel > 0) {
	r1.push_back (true);
	r2.push_back (false);
	--nDel;
      }
    case ProbModel::Insert:
      r1.push_back (c1);
      r2.push_back (c2);
      break;
    case ProbModel::Delete:
      ++nDel;
    default:
      break;
    }
  }
  while (nDel > 0) {
    r1.push_back (true);
    r2.push_back (false);
    --nDel;
  }
  Assert (alignPathResiduesInRow(r1) == alignPathResiduesInRow(row1)
	  && alignPathResiduesInRow(r2) == alignPathResiduesInRow(row2), "Rows don't match");
  return p;
}

AlignPath TreeAlignFuncs::triplePath (const AlignPath& path, TreeNodeIndex lChild, TreeNodeIndex rChild, TreeNodeIndex parent) {
  AlignPath p;
  size_t nLeftIns = 0;
  const AlignColIndex cols = alignPathColumns (path);
  const AlignRowPath& lRow = path.at (lChild);
  const AlignRowPath& rRow = path.at (rChild);
  const AlignRowPath& pRow = path.at (parent);
  AlignRowPath& lr = p[lChild];
  AlignRowPath& rr = p[rChild];
  AlignRowPath& pr = p[parent];
  for (AlignColIndex col = 0; col < cols; ++col) {
    const bool lc = lRow[col];
    const bool rc = rRow[col];
    const bool pc = pRow[col];
    if (!(lc || rc || pc))
      continue;
    const auto state = Sampler::SiblingMatrix::getState (Sampler::SiblingMatrix::IMM, lc, rc, pc);
    switch (state) {
    case Sampler::SiblingMatrix::IMM:
    case Sampler::SiblingMatrix::IMD:
    case Sampler::SiblingMatrix::IDM:
    case Sampler::SiblingMatrix::IDD:
      while (nLeftIns > 0) {
	lr.push_back (true);
	rr.push_back (false);
	pr.push_back (false);
	--nLeftIns;
      }
    case Sampler::SiblingMatrix::IMI:
      lr.push_back (lc);
      rr.push_back (rc);
      pr.push_back (pc);
      break;
    case Sampler::SiblingMatrix::IIW:
      ++nLeftIns;
      break;
    default:
      Abort ("bad state: %d  (l,r,p)=(%d,%d,%d)", (int) state, lc, rc, pc);
      break;
    }
  }
  while (nLeftIns > 0) {
    lr.push_back (true);
    rr.push_back (false);
    pr.push_back (false);
    --nLeftIns;
  }
  Assert (alignPathResiduesInRow(lr) == alignPathResiduesInRow(lRow)
	  && alignPathResiduesInRow(rr) == alignPathResiduesInRow(rRow)
	  && alignPathResiduesInRow(pr) == alignPathResiduesInRow(pRow), "Rows don't match");
  return p;
}

AlignPath TreeAlignFuncs::branchPath (const AlignPath& path, const Tree& tree, TreeNodeIndex node) {
  const TreeNodeIndex parent = tree.parentNode(node);
  Assert (parent >= 0, "Parent node not found");
  return pairPath (path, parent, node);
}

bool TreeAlignFuncs::subpathUngapped (const AlignPath& path, const vguard<TreeNodeIndex>& rows) {
  const AlignColIndex cols = alignPathColumns (path);
  for (AlignColIndex col = 0; col < cols; ++col) {
    size_t n = 0;
    for (auto row : rows)
      if (path.at(row)[col])
	++n;
    if (n != 0 && n != rows.size())
      return false;
  }
  return true;
}

set<TreeNodeIndex> TreeAlignFuncs::allNodes (const Tree& tree) {
  const auto nvec = tree.preorderSort();
  return set<TreeNodeIndex> (nvec.begin(), nvec.end());
}

set<TreeNodeIndex> TreeAlignFuncs::allExceptNodeAndAncestors (const Tree& tree, TreeNodeIndex node) {
  set<TreeNodeIndex> n = allNodes (tree), a = tree.nodeAndAncestors (node);
  for (auto anc: a)
    n.erase (anc);
  return n;
}

set<TreeNodeIndex> TreeAlignFuncs::nodeAndAncestors (const Tree& tree, TreeNodeIndex node) {
  return tree.nodeAndAncestors (node);
}

set<TreeNodeIndex> TreeAlignFuncs::nodesAndAncestors (const Tree& tree, TreeNodeIndex node1, TreeNodeIndex node2) {
  set<TreeNodeIndex> a = tree.nodeAndAncestors (node1), n2 = tree.nodeAndAncestors (node2);
  a.insert (n2.begin(), n2.end());
  return a;
}

string TreeAlignFuncs::branchConditionalDump (const RateModel& model, const Tree& tree, const vguard<FastSeq>& gapped, TreeNodeIndex parent, TreeNodeIndex node) {
  ostringstream out;

  const ProbModel probModel (model, max (Tree::minBranchLength, tree.branchLength(parent,node)));
  const LogProbModel logProbModel (probModel);
  const auto& submat (logProbModel.logSubProb);

  map<TreeNodeIndex,TreeNodeIndex> exclude;
  exclude[node] = parent;
  exclude[parent] = node;

  const set<TreeNodeIndex> fillUpNodes = allExceptNodeAndAncestors (tree, parent);
  const set<TreeNodeIndex> fillDownNodes = nodeAndAncestors (tree, parent);

  const auto pwms = getConditionalPWMs (model, tree, gapped, exclude, fillUpNodes, fillDownNodes, false);
  const PosWeightMatrix& pSeq = pwms.at (parent);
  const PosWeightMatrix& nSeq = pwms.at (node);
  const PosWeightMatrix nPre = preMultiply (nSeq, submat);

  AlignColSumProduct colSumProdBranch (model, tree, gapped);
  AlignColSumProduct colSumProdFull (model, tree, gapped);

  colSumProdBranch.preorder = vguard<TreeNodeIndex> (fillDownNodes.rbegin(), fillDownNodes.rend());
  colSumProdBranch.postorder = vguard<TreeNodeIndex> (fillUpNodes.begin(), fillUpNodes.end());

  size_t col = 0, pCol = 0, nCol = 0;
  while (!colSumProdBranch.alignmentDone()) {
    colSumProdBranch.fillUp();
    colSumProdBranch.fillDown();
    colSumProdFull.fillUp();
    colSumProdFull.fillDown();

    const bool pGap = Alignment::isGap(gapped[parent].seq[col]);
    const bool nGap = Alignment::isGap(gapped[node].seq[col]);

    const LogProb cll = colSumProdFull.columnLogLikelihood();
    const LogProb pll = pGap ? -numeric_limits<double>::infinity() : colSumProdFull.computeColumnLogLikelihoodAt (parent);
    const LogProb nll = nGap ? -numeric_limits<double>::infinity() : colSumProdFull.computeColumnLogLikelihoodAt (node);

    LogProb lp = -numeric_limits<double>::infinity();
    LogProb lpc = -numeric_limits<double>::infinity();
    if (!pGap && !nGap) {
      const auto pCond = colSumProdBranch.logNodeExcludedPostProb (parent, node, false);
      const auto nCond = colSumProdBranch.logNodeExcludedPostProb (node, parent, false);

      for (int cpt = 0; cpt < model.components(); ++cpt)
	for (AlphTok i = 0; i < model.alphabetSize(); ++i) {
	  log_accum_exp (lp, pSeq[pCol][cpt][i] + nPre[nCol][cpt][i]);
	  for (AlphTok j = 0; j < model.alphabetSize(); ++j)
	    log_accum_exp (lpc, pCond[cpt][i] + submat[cpt][i][j] + nCond[cpt][j]);
	}
    }

    out << "Column #" << col << " (parent #" << pCol << ", node #" << nCol << ") "
	<< "log-likelihood: " << cll << "(root) " << pll << "(parent) " << nll << "(node) "
	<< lp << "(pwm) " << lpc << "(branch)"
	<< endl;

    colSumProdBranch.nextColumn();
    colSumProdFull.nextColumn();

    if (!pGap)
      ++pCol;
    if (!nGap)
      ++nCol;
    ++col;
  }

  return out.str();
}

map<TreeNodeIndex,TreeAlignFuncs::PosWeightMatrix> TreeAlignFuncs::getConditionalPWMs (const RateModel& model, const Tree& tree, const vguard<FastSeq>& gapped, const map<TreeNodeIndex,TreeNodeIndex>& exclude, const set<TreeNodeIndex>& fillUpNodes, const set<TreeNodeIndex>& fillDownNodes, bool normalize) {
  map<TreeNodeIndex,PosWeightMatrix> pwms;
  AlignColSumProduct colSumProd (model, tree, gapped);
  colSumProd.preorder = vguard<TreeNodeIndex> (fillDownNodes.rbegin(), fillDownNodes.rend());
  colSumProd.postorder = vguard<TreeNodeIndex> (fillUpNodes.begin(), fillUpNodes.end());
  while (!colSumProd.alignmentDone()) {
    colSumProd.fillUp();
    colSumProd.fillDown();
    for (const auto& node_exclude : exclude)
      if (!colSumProd.isGap (node_exclude.first))
	pwms[node_exclude.first].push_back (colSumProd.logNodeExcludedPostProb (node_exclude.first, node_exclude.second, normalize));
    colSumProd.nextColumn();
  }
  return pwms;
}

LogProb TreeAlignFuncs::rootLogLikelihood (const RateModel& model, const History& history) {
  size_t rootLen = 0;
  for (auto c: history.gapped[history.tree.root()].seq)
    if (!Alignment::isGap(c))
      ++rootLen;
  const LogProb rootExt = rootExtProb(model);
  const LogProb lpRoot = log(1. - rootExt) + log(rootExt) * rootLen;
  return lpRoot;
}

LogProb TreeAlignFuncs::indelLogLikelihood (const RateModel& model, const History& history) {
  const Alignment align (history.gapped);
  LogProb lpGaps = 0;
  for (TreeNodeIndex node = 0; node < history.tree.root(); ++node) {
    const TreeNodeIndex parent = history.tree.parentNode (node);
    const ProbModel probModel (model, history.tree.branchLength (node));
    const AlignPath path = pairPath (align.path, parent, node);
    lpGaps += logBranchPathLikelihood (probModel, path, parent, node);
  }
  return lpGaps;
}

LogProb TreeAlignFuncs::substLogLikelihood (const RateModel& model, const History& history) {
  AlignColSumProduct colSumProd (model, history.tree, history.gapped);
  LogProb lpSub = 0;
  vguard<LogProb> colSub;
  while (!colSumProd.alignmentDone()) {
    colSumProd.fillUp();
    const LogProb cll = colSumProd.columnLogLikelihood();
    lpSub += cll;
    colSub.push_back (cll);
    colSumProd.nextColumn();
  }
  LogThisAt(9,"Column substitution log-likelihoods: (" << to_string_join(colSub) << ")" << endl);
  return lpSub;
}

LogProb TreeAlignFuncs::logLikelihood (const SimpleTreePrior& treePrior, const RateModel& model, const History& history, const char* suffix) {
  const LogProb lpTree = treePrior.treeLogLikelihood (history.tree);
  const LogProb lpRoot = rootLogLikelihood (model, history);
  const LogProb lpGaps = indelLogLikelihood (model, history);
  const LogProb lpSub = substLogLikelihood (model, history);
  const LogProb lp = lpTree + lpRoot + lpGaps + lpSub;
  LogThisAt(6,"log(L" << suffix << ") = " << setw(10) << lpTree << " (tree) + " << setw(10) << lpRoot << " (root) + " << setw(10) << lpGaps << " (indels) + " << setw(10) << lpSub << " (substitutions) = " << lp << endl);
  return lp;
}

LogProb TreeAlignFuncs::logLikelihood (const RateModel& model, const History& history, const char* suffix) {
  const Alignment align (history.gapped);
  const LogProb lpRoot = rootLogLikelihood (model, history);
  const LogProb lpGaps = indelLogLikelihood (model, history);
  const LogProb lpSub = substLogLikelihood (model, history);
  const LogProb lp = lpRoot + lpGaps + lpSub;
  LogThisAt(6,"log(L" << suffix << ") = " << setw(10) << lpRoot << " (root) + " << setw(10) << lpGaps << " (indels) + " << setw(10) << lpSub << " (substitutions) = " << lp << endl);
  return lp;
}

LogProb TreeAlignFuncs::logLikelihood (const SimpleTreePrior& treePrior, const RateModel& model, const Tree& tree, const vguard<FastSeq>& gapped, const char* suffix) {
  const History history (gapped, tree);
  return logLikelihood (treePrior, model, history, suffix);
}

LogProb TreeAlignFuncs::logLikelihood (const RateModel& model, const Tree& tree, const vguard<FastSeq>& gapped, const char* suffix) {
  const History history (gapped, tree);
  return logLikelihood (model, history, suffix);
}

LogProb TreeAlignFuncs::logBranchPathLikelihood (const ProbModel& probModel, const AlignPath& path, TreeNodeIndex parent, TreeNodeIndex child) {
  const AlignColIndex cols = alignPathColumns (path);
  ProbModel::State state = ProbModel::Start;
  LogProb lp = 0;
  for (AlignColIndex col = 0; col < cols; ++col) {
    const ProbModel::State nextState = ProbModel::getState (path.at(parent)[col], path.at(child)[col]);
    lp += log (probModel.transProb (state, nextState));
    state = nextState;
  }
  lp += log (probModel.transProb (state, ProbModel::End));
  return lp;
}

TreeAlignFuncs::PosWeightMatrix TreeAlignFuncs::preMultiply (const PosWeightMatrix& child, const vguard<LogProbModel::LogProbMatrix>& submat) {
  PosWeightMatrix pwm (child.size(), vguard<vguard<LogProb> > (submat.size(), vguard<LogProb> (submat.front().size(), -numeric_limits<double>::infinity())));
  size_t n = 0;
  for (const auto& lpp : child) {
    auto& pre = pwm[n++];
    for (int cpt = 0; cpt < (int) submat.size(); ++cpt)
      for (AlphTok i = 0; i < pre[cpt].size(); ++i)
	for (AlphTok j = 0; j < lpp[cpt].size(); ++j)
	  log_accum_exp (pre[cpt][i], submat[cpt][i][j] + lpp[cpt][j]);
  }
  return pwm;
}

vguard<LogProb> TreeAlignFuncs::calcInsProbs (const PosWeightMatrix& child, const vguard<LogProbModel::LogProbVector>& insvec, const vguard<LogProb>& logCptWeight) {
  vguard<LogProb> ins;
  ins.reserve (child.size());
  for (const auto& lpp : child) {
    LogProb lp = -numeric_limits<double>::infinity();
    for (int cpt = 0; cpt < (int) insvec.size(); ++cpt)
      for (AlphTok i = 0; i < lpp[cpt].size(); ++i)
	log_accum_exp (lp, logCptWeight[cpt] + insvec[cpt][i] + lpp[cpt][i]);
    ins.push_back (lp);
  }
  return ins;
}

Sampler::Move::Move (Type type, const History& history, LogProb oldLogLikelihood, const string& samplerName)
  : type (type),
    nullified (false),
    newLogLikelihood (0),
    oldLogLikelihood (oldLogLikelihood),
    logForwardProposal (0),
    logReverseProposal (0),
    logJacobian (0),
    logAcceptProb (-numeric_limits<double>::infinity()),
    oldHistory (history),
    samplerName (samplerName)
{ }

void Sampler::Move::initNewHistory (const Tree& tree, const vguard<FastSeq>& ungapped, const AlignPath& path) {
  const Alignment newAlign (ungapped, path);
  initNewHistory (tree, newAlign.gapped());
}

void Sampler::Move::initNewHistory (const Tree& tree, const vguard<FastSeq>& gapped) {
  newHistory.tree = tree;
  newHistory.gapped = gapped;
}

void Sampler::Move::initRatio (const Sampler& sampler) {
  newLogLikelihood = sampler.logLikelihood (newHistory, "_new");

  const LogProb logOddsRatio = newLogLikelihood - oldLogLikelihood;
  const LogProb logHastingsRatio = logReverseProposal - logForwardProposal + logJacobian;
  logAcceptProb = logOddsRatio + logHastingsRatio;

  LogThisAt(5,"log(L_new/L_old) = " << logOddsRatio << ", log(Q_rev/Q_fwd) = " << logHastingsRatio << ", log(P_accept) = " << logAcceptProb << endl);
}

void Sampler::Move::nullify (const char* reason) {
  newHistory = oldHistory;
  newLogLikelihood = oldLogLikelihood;
  logAcceptProb = logJacobian = logForwardProposal = logReverseProposal = 0;
  nullified = true;
  comment = string("(") + reason + ")";
}

const char* Sampler::Move::typeName (Type t) {
  switch (t) {
  case BranchAlign:
    return "Branch alignment";
  case NodeAlign:
    return "Node alignment";
  case PruneAndRegraft:
    return "Prune-and-regraft";
  case NodeHeight:
    return "Node height";
  case Rescale:
    return "Rescale";
  default:
    break;
  }
  return "Unknown";
}

size_t Sampler::Move::typeNameWidth() {
  return 20;
}

bool Sampler::Move::accept (random_engine& generator) const {
  bool a;
  if (nullified)
    a = false;
  else if (logAcceptProb >= 0)
    a = true;
  else {
    bernoulli_distribution distribution (exp (logAcceptProb));
    a = distribution (generator);
  }
  LogThisAt(3,setw(30) << left << samplerName << right << " "
	    << setw(Move::typeNameWidth()) << typeName(type) << " move "
	    << (nullified ? "bypassed" : (a ? "ACCEPTED" : "rejected"))
	    << " with log(P_accept) = " << setw(10) << logAcceptProb
	    << " " << comment << endl);
  return a;
}

Sampler::BranchAlignMove::BranchAlignMove (const History& history, LogProb oldLogLikelihood, const Sampler& sampler, random_engine& generator)
  : Move (BranchAlign, history, oldLogLikelihood, sampler.name)
{
  node = Sampler::randomChildNode (history.tree, generator);
  parent = history.tree.parentNode (node);

  LogThisAt(4,"Proposing branch realignment move between...\n   node #" << node << ": " << history.tree.seqName(node) << "\n parent #" << parent << ": " << history.tree.seqName(parent) << endl);

  const TreeBranchLength dist = history.tree.branchLength(parent,node);
  
  const TreeNodeIndex parentClosestLeaf = history.tree.closestLeaf (parent, node);
  const TreeNodeIndex nodeClosestLeaf = history.tree.closestLeaf (node, parent);

  const Alignment oldAlign (history.gapped);
  const AlignPath oldBranchPath = Sampler::branchPath (oldAlign.path, history.tree, node);
  const GuideAlignmentEnvelope newBranchEnv = sampler.makeGuide (history.tree, parentClosestLeaf, nodeClosestLeaf, oldBranchPath, parent, node);

  const AlignPath pCladePath = Sampler::cladePath (oldAlign.path, history.tree, parent, node);
  const AlignPath nCladePath = Sampler::cladePath (oldAlign.path, history.tree, node, parent);
  
  const vguard<SeqIdx> parentEnvPos = sampler.guideSeqPos (oldAlign.path, parent, parentClosestLeaf);
  const vguard<SeqIdx> nodeEnvPos = sampler.guideSeqPos (oldAlign.path, node, nodeClosestLeaf);

  map<TreeNodeIndex,TreeNodeIndex> exclude;
  exclude[node] = parent;
  exclude[parent] = node;

  const auto pwms = sampler.getConditionalPWMs (history.tree, history.gapped, exclude, sampler.allExceptNodeAndAncestors(history.tree,parent), sampler.nodeAndAncestors(history.tree,parent));
  const PosWeightMatrix& pSeq = pwms.at (parent);
  const PosWeightMatrix& nSeq = pwms.at (node);

  const BranchMatrix newBranchMatrix (sampler.model, pSeq, nSeq, dist, newBranchEnv, parentEnvPos, nodeEnvPos, parent, node);
  const AlignPath newBranchPath = newBranchMatrix.sample (generator);

  LogThisAt(6,"Proposed (parent:node) alignment:" << endl << alignPathString(newBranchPath));
  const LogProb logPostNewBranchPath = newBranchMatrix.logPostProb (newBranchPath);

  LogThisAt(6,"Previous (parent:node) alignment:" << endl << alignPathString(oldBranchPath));
  const GuideAlignmentEnvelope oldBranchEnv = sampler.makeGuide (history.tree, parentClosestLeaf, nodeClosestLeaf, newBranchPath, parent, node);
  const BranchMatrix* oldBranchMatrix =
    sampler.useFixedGuide
    ? &newBranchMatrix
    : new BranchMatrix (sampler.model, pSeq, nSeq, dist, oldBranchEnv, parentEnvPos, nodeEnvPos, parent, node);
  const LogProb logPostOldBranchPath = oldBranchMatrix->logPostProb (oldBranchPath);
  if (!sampler.useFixedGuide)
    delete oldBranchMatrix;
  
  if (oldBranchPath == newBranchPath) {
    LogThisAt(6,"Alignments are identical; abandoning move" << endl);
    nullify("no change");
    return;
  }

  const vguard<AlignPath> mergeComponents = { pCladePath, newBranchPath, nCladePath };
  const AlignPath newPath = alignPathMerge (mergeComponents);

  logForwardProposal = logPostNewBranchPath;
  logReverseProposal = logPostOldBranchPath;

  initNewHistory (oldHistory.tree, oldAlign.ungapped, newPath);
  initRatio (sampler);
}

Sampler::NodeAlignMove::NodeAlignMove (const History& history, LogProb oldLogLikelihood, const Sampler& sampler, random_engine& generator)
  : Move (NodeAlign, history, oldLogLikelihood, sampler.name)
{
  node = Sampler::randomInternalNode (history.tree, generator);

  Assert (history.tree.nChildren(node) == 2, "Non-binary tree");
  leftChild = history.tree.getChild (node, 0);
  rightChild = history.tree.getChild (node, 1);

  parent = history.tree.parentNode (node);

  LogThisAt(4,"Proposing node realignment move between...\n        node #" << node << ": " << history.tree.seqName(node) << "\n  left-child #" << leftChild << ": " << history.tree.seqName(leftChild) << "\n right-child #" << rightChild << ": " << history.tree.seqName(rightChild) << (parent >= 0 ? (string("\n      parent #") + to_string(parent) + ": " + history.tree.seqName(parent)) : string()) << endl);

  const TreeBranchLength lDist = history.tree.branchLength(node,leftChild);
  const TreeBranchLength rDist = history.tree.branchLength(node,rightChild);
  
  const TreeNodeIndex leftChildClosestLeaf = history.tree.closestLeaf (leftChild, node);
  const TreeNodeIndex rightChildClosestLeaf = history.tree.closestLeaf (rightChild, node);

  const Alignment oldAlign (history.gapped);
  const AlignPath oldSiblingPath = Sampler::triplePath (oldAlign.path, leftChild, rightChild, node);

  const AlignPath lCladePath = Sampler::cladePath (oldAlign.path, history.tree, leftChild, node);
  const AlignPath rCladePath = Sampler::cladePath (oldAlign.path, history.tree, rightChild, node);

  const vguard<SeqIdx> leftChildEnvPos = sampler.guideSeqPos (oldAlign.path, leftChild, leftChildClosestLeaf);
  const vguard<SeqIdx> rightChildEnvPos = sampler.guideSeqPos (oldAlign.path, rightChild, rightChildClosestLeaf);

  const GuideAlignmentEnvelope newSiblingEnv = sampler.makeGuide (history.tree, leftChildClosestLeaf, rightChildClosestLeaf, oldSiblingPath, leftChild, rightChild);

  map<TreeNodeIndex,TreeNodeIndex> exclude;
  exclude[leftChild] = node;
  exclude[rightChild] = node;
  if (parent >= 0) {
    exclude[node] = parent;
    exclude[parent] = node;
  }
  const auto pwms = sampler.getConditionalPWMs (history.tree, history.gapped, exclude, sampler.allExceptNodeAndAncestors(history.tree,parent>=0?parent:node), parent >= 0 ? sampler.nodeAndAncestors(history.tree,parent) : set<TreeNodeIndex>());
  const PosWeightMatrix& lSeq = pwms.at (leftChild);
  const PosWeightMatrix& rSeq = pwms.at (rightChild);

  const SiblingMatrix newSibMatrix (sampler.model, lSeq, rSeq, lDist, rDist, newSiblingEnv, leftChildEnvPos, rightChildEnvPos, leftChild, rightChild, node);
  const AlignPath newSiblingPath = newSibMatrix.sample (generator);
  LogThisAt(6,"Proposed (node:left:right) alignment:" << endl << alignPathString(newSiblingPath));
  const LogProb logPostNewSiblingPath = newSibMatrix.logPostProb (newSiblingPath);

  LogThisAt(6,"Previous (node:left:right) alignment:" << endl << alignPathString(oldSiblingPath));
  const GuideAlignmentEnvelope oldSiblingEnv = sampler.makeGuide (history.tree, leftChildClosestLeaf, rightChildClosestLeaf, newSiblingPath, leftChild, rightChild);
  const SiblingMatrix* oldSibMatrix =
    sampler.useFixedGuide
    ? &newSibMatrix
    : new SiblingMatrix (sampler.model, lSeq, rSeq, lDist, rDist, oldSiblingEnv, leftChildEnvPos, rightChildEnvPos, leftChild, rightChild, node);
  const LogProb logPostOldSiblingPath = oldSibMatrix->logPostProb (oldSiblingPath);

  logForwardProposal = logPostNewSiblingPath;
  logReverseProposal = logPostOldSiblingPath;

  vguard<AlignPath> mergeComponents = { lCladePath, rCladePath, newSiblingPath };
  AlignPath newPath = alignPathMerge (mergeComponents);

  const PosWeightMatrix newNodeSeq = newSibMatrix.parentSeq (newSiblingPath);
  const PosWeightMatrix& oldNodeSeq = oldSibMatrix->parentSeq (oldSiblingPath);

  const vguard<FastSeq>& oldUngapped = oldAlign.ungapped;
  vguard<FastSeq> newUngapped = oldUngapped;
  if (sampler.sampleAncestralSeqs) {
    newUngapped[node].seq = sampler.sampleSeq (newNodeSeq, generator);
    logForwardProposal += sampler.logSeqPostProb (newUngapped[node].seq, newNodeSeq);
    logReverseProposal += sampler.logSeqPostProb (oldUngapped[node].seq, oldNodeSeq);
  } else
    newUngapped[node].seq = string (alignPathResiduesInRow (newSiblingPath.at (node)), Alignment::wildcardChar);

  if (parent >= 0) {  // don't attempt to align to parent if node is root
    const PosWeightMatrix& pSeq = pwms.at (parent);
    const TreeBranchLength pDist = history.tree.branchLength(parent,node);

    const TreeNodeIndex nodeClosestLeaf = history.tree.closestLeaf (node, parent);
    const TreeNodeIndex parentClosestLeaf = history.tree.closestLeaf (parent, node);
    const TreeNodeIndex nodeClosestChild = lDist < rDist ? leftChild : rightChild;
    
    const GuideAlignmentEnvelope newBranchEnv = sampler.makeGuide (history.tree, parentClosestLeaf, nodeClosestLeaf, oldAlign.path, parent, nodeClosestChild);
    const AlignPath& nodeSubtreePath = newPath;
    const vguard<SeqIdx> newNodeEnvPos = sampler.guideSeqPos (nodeSubtreePath, node, nodeClosestChild, nodeClosestLeaf);
    const vguard<SeqIdx> oldNodeEnvPos = sampler.guideSeqPos (oldAlign.path, node, nodeClosestChild, nodeClosestLeaf);

    const AlignPath pCladePath = Sampler::cladePath (oldAlign.path, history.tree, parent, node);
    const vguard<SeqIdx> parentEnvPos = sampler.guideSeqPos (oldAlign.path, parent, parentClosestLeaf);

    const BranchMatrix newBranchMatrix (sampler.model, pSeq, newNodeSeq, pDist, newBranchEnv, parentEnvPos, newNodeEnvPos, parent, node);

    const AlignPath newBranchPath = newBranchMatrix.sample (generator);
    LogThisAt(6,"Proposed (parent:node) alignment:" << endl << alignPathString(newBranchPath));
    const LogProb logPostNewBranchPath = newBranchMatrix.logPostProb (newBranchPath);

    mergeComponents.push_back (pCladePath);
    mergeComponents.push_back (newBranchPath);
    newPath = alignPathMerge (mergeComponents);

    const GuideAlignmentEnvelope oldBranchEnv = sampler.makeGuide (history.tree, parentClosestLeaf, nodeClosestLeaf, newPath, parent, nodeClosestChild);
    const BranchMatrix oldBranchMatrix (sampler.model, pSeq, oldNodeSeq, pDist, oldBranchEnv, parentEnvPos, oldNodeEnvPos, parent, node);
    const AlignPath oldBranchPath = Sampler::branchPath (oldAlign.path, history.tree, node);
    LogThisAt(6,"Previous (parent:node) alignment:" << endl << alignPathString(oldBranchPath));
    const LogProb logPostOldBranchPath = oldBranchMatrix.logPostProb (oldBranchPath);

    logForwardProposal += logPostNewBranchPath;
    logReverseProposal += logPostOldBranchPath;

    LogThisAt(6,"log(Q_fwd) = " << setw(10) << logPostNewSiblingPath << " (node:leftChild:rightChild) + " << setw(10) << logPostNewBranchPath << " (parent:node) = " << logForwardProposal << endl);
    LogThisAt(6,"log(Q_rev) = " << setw(10) << logPostOldSiblingPath << " (node:leftChild:rightChild) + " << setw(10) << logPostOldBranchPath << " (parent:node) = " << logReverseProposal << endl);
  }

  if (!sampler.useFixedGuide)
    delete oldSibMatrix;

  if (newPath == oldAlign.path && (!sampler.sampleAncestralSeqs || newUngapped[node].seq == oldUngapped[node].seq)) {
    LogThisAt(6,"Alignments are identical; abandoning move" << endl);
    nullify("no change");
    return;
  }
  
  initNewHistory (oldHistory.tree, newUngapped, newPath);
  initRatio (sampler);
}

Sampler::PruneAndRegraftMove::PruneAndRegraftMove (const History& history, LogProb oldLogLikelihood, const Sampler& sampler, random_engine& generator)
  : Move (PruneAndRegraft, history, oldLogLikelihood, sampler.name)
{
  const vguard<TreeBranchLength>& distanceFromRoot = history.tree.distanceFromRoot();

  node = Sampler::randomGrandchildNode (history.tree, generator);

  const vector<TreeNodeIndex> contemps = contemporaneousNodes (history.tree, distanceFromRoot, node);
  if (contemps.empty()) {
    nullify("nowhere to regraft");
    return;
  }
  const vguard<double> contempWeights = sampler.nodeListWeights (contemps.size());
  const size_t contempIndex = random_index (contempWeights, generator);
  const TreeNodeIndex newSibling = contemps[contempIndex];

  parent = history.tree.parentNode (node);
  Assert (parent >= 0, "Parent node not found");

  oldGrandparent = history.tree.parentNode (parent);
  Assert (oldGrandparent >= 0, "Grandparent node not found");

  newGrandparent = history.tree.parentNode (newSibling);
  Assert (newGrandparent >= 0, "Grandparent node not found");

  oldSibling = history.tree.getSibling (node);
  
  LogThisAt(4,"Proposing prune-and-regraft move at...\n            node #" << node << ": " << history.tree.seqName(node) << "\n          parent #" << parent << ": " << history.tree.seqName(parent) << "\n     old sibling #" << oldSibling << ": " << history.tree.seqName(oldSibling) << "\n old grandparent #" << oldGrandparent << ": " << history.tree.seqName(oldGrandparent) << "\n     new sibling #" << newSibling << ": " << history.tree.seqName(newSibling) << "\n new grandparent #" << newGrandparent << ": " << history.tree.seqName(newGrandparent) << endl);

  const Tree& oldTree = history.tree;
  const Alignment oldAlign (history.gapped);
  
  const TreeBranchLength oldGranParentDist = oldTree.branchLength(oldGrandparent,parent);
  const TreeBranchLength parentNodeDist = oldTree.branchLength(parent,node);
  const TreeBranchLength parentOldSibDist = oldTree.branchLength(parent,oldSibling);

  const TreeBranchLength parentNewSibDist = distanceFromRoot[newSibling] - distanceFromRoot[parent];
  const TreeBranchLength newGranParentDist = distanceFromRoot[parent] - distanceFromRoot[newGrandparent];

  Tree newTree = history.tree;
  newTree.setParent (oldSibling, oldGrandparent, oldGranParentDist + parentOldSibDist);
  newTree.setParent (newSibling, parent, parentNewSibDist);
  newTree.setParent (parent, newGrandparent, newGranParentDist);

  const vector<TreeNodeIndex> revContemps = contemporaneousNodes (newTree, newTree.distanceFromRoot(), node);
  const vguard<double> revContempWeights = sampler.nodeListWeights (revContemps.size());
  const size_t revContempIndex = find (revContemps.begin(), revContemps.end(), oldSibling) - revContemps.begin();
  if (revContempIndex >= revContemps.size()) {
    Warn ("Couldn't invert move");
    nullify("couldn't invert move");
    return;
  }

  const LogProb logFwdSibSelect = log(contempWeights[contempIndex]);
  const LogProb logRevSibSelect = log(revContempWeights[revContempIndex]);

  // optimize special case that (oldSibling,parent,oldGrandparent,newGrandparent,newSibling) form a sub-alignment with no gaps
  const vguard<TreeNodeIndex> subpathNodes = { oldSibling, parent, oldGrandparent, newGrandparent, newSibling };
  if (Sampler::subpathUngapped (oldAlign.path, subpathNodes)) {
    initNewHistory (newTree, history.gapped);

    logForwardProposal = logFwdSibSelect;
    logReverseProposal = logRevSibSelect;

    comment = "(alignment unchanged)";
    
  } else {
    // general case: we need to realign
    const AlignPath oldSibCladePath = Sampler::cladePath (oldAlign.path, oldTree, oldSibling, parent);
    const AlignPath nodeCladePath = Sampler::cladePath (oldAlign.path, oldTree, node, parent);
    const AlignPath newSibCladePath = Sampler::cladePath (oldAlign.path, oldTree, newSibling, newGrandparent);
    const AlignPath oldGranCladePath = Sampler::cladePath (oldAlign.path, oldTree, oldGrandparent, parent, newSibling);

    const AlignPath oldSiblingPath = Sampler::triplePath (oldAlign.path, node, oldSibling, parent);
    const AlignPath oldBranchPath = Sampler::branchPath (oldAlign.path, oldTree, parent);

    const AlignPath oldGranSibPath = Sampler::pairPath (oldAlign.path, oldGrandparent, oldSibling);
    
    const TreeNodeIndex nodeClosestLeaf = oldTree.closestLeaf (node, parent);
    const TreeNodeIndex oldSibClosestLeaf = oldTree.closestLeaf (oldSibling, parent);
    const TreeNodeIndex oldGranClosestLeaf = oldTree.closestLeaf (oldGrandparent, parent);
    const TreeNodeIndex newSibClosestLeaf = newTree.closestLeaf (newSibling, parent);
    const TreeNodeIndex newGranClosestLeaf = newTree.closestLeaf (newGrandparent, parent);
    const TreeNodeIndex oldParentClosestLeaf = oldTree.closestLeaf (parent, oldGrandparent);
    const TreeNodeIndex newParentClosestLeaf = newTree.closestLeaf (parent, newGrandparent);

    const TreeNodeIndex oldParentClosestChild = parentNodeDist < parentOldSibDist ? node : oldSibling;
    const TreeNodeIndex newParentClosestChild = parentNodeDist < parentNewSibDist ? node : newSibling;

    const vguard<SeqIdx> nodeEnvPos = sampler.guideSeqPos (oldAlign.path, node, nodeClosestLeaf);
    const vguard<SeqIdx> oldSibEnvPos = sampler.guideSeqPos (oldAlign.path, oldSibling, oldSibClosestLeaf);
    const vguard<SeqIdx> oldGranEnvPos = sampler.guideSeqPos (oldAlign.path, oldGrandparent, oldGranClosestLeaf);
    const vguard<SeqIdx> newSibEnvPos = sampler.guideSeqPos (oldAlign.path, newSibling, newSibClosestLeaf);
    const vguard<SeqIdx> newGranEnvPos = sampler.guideSeqPos (oldAlign.path, newGrandparent, newGranClosestLeaf);

    const GuideAlignmentEnvelope newSibEnv = sampler.makeGuide (history.tree, nodeClosestLeaf, newSibClosestLeaf, oldAlign.path, node, newSibling);

    map<TreeNodeIndex,TreeNodeIndex> exclude;
    exclude[node] = -1;
    exclude[oldSibling] = parent;
    exclude[oldGrandparent] = parent;
    exclude[newSibling] = newGrandparent;
    exclude[newGrandparent] = newSibling;

    Tree detachedTree = oldTree;
    detachedTree.detach (node);
    const auto pwms = sampler.getConditionalPWMs (detachedTree, history.gapped, exclude, sampler.allNodes(history.tree), sampler.nodesAndAncestors(history.tree,oldGrandparent,newGrandparent));

    const PosWeightMatrix& nodeSeq = pwms.at (node);
    const PosWeightMatrix& oldSibSeq = pwms.at (oldSibling);
    const PosWeightMatrix& oldGranSeq = pwms.at (oldGrandparent);
    const PosWeightMatrix& newSibSeq = pwms.at (newSibling);
    const PosWeightMatrix& newGranSeq = pwms.at (newGrandparent);
    
    const SiblingMatrix newSibMatrix (sampler.model, nodeSeq, newSibSeq, parentNodeDist, parentNewSibDist, newSibEnv, nodeEnvPos, newSibEnvPos, node, newSibling, parent);
    const AlignPath newSiblingPath = newSibMatrix.sample (generator);
    LogThisAt(6,"Proposed (parent:node:newSibling) alignment:" << endl << alignPathString(newSiblingPath));
    const LogProb logPostNewSiblingPath = newSibMatrix.logPostProb (newSiblingPath);

    LogThisAt(6,"Previous (parent:node:oldSibling) alignment:" << endl << alignPathString(oldSiblingPath));

    vguard<AlignPath> mergeComponents = { nodeCladePath, newSibCladePath, newSiblingPath };
    const AlignPath newParentSubtreePath = alignPathMerge (mergeComponents);

    const GuideAlignmentEnvelope newBranchEnv = sampler.makeGuide (history.tree, newGranClosestLeaf, newParentClosestLeaf, oldAlign.path, newGrandparent, newParentClosestChild);

    const vguard<SeqIdx> newParentEnvPos = sampler.guideSeqPos (newParentSubtreePath, parent, newParentClosestChild, newParentClosestLeaf);
    const vguard<SeqIdx> oldParentEnvPos = sampler.guideSeqPos (oldAlign.path, parent, oldParentClosestChild, oldParentClosestLeaf);

    const PosWeightMatrix newParentSeq = newSibMatrix.parentSeq (newSiblingPath);
    const BranchMatrix newBranchMatrix (sampler.model, newGranSeq, newParentSeq, newGranParentDist, newBranchEnv, newGranEnvPos, newParentEnvPos, newGrandparent, parent);

    const AlignPath newBranchPath = newBranchMatrix.sample (generator);
    LogThisAt(6,"Proposed (newGrandparent:parent) alignment:" << endl << alignPathString(newBranchPath));
    const LogProb logPostNewBranchPath = newBranchMatrix.logPostProb (newBranchPath);

    LogThisAt(6,"Previous (oldGrandparent:parent) alignment:" << endl << alignPathString(oldBranchPath));

    mergeComponents.push_back (oldSibCladePath);
    mergeComponents.push_back (oldGranSibPath);
    mergeComponents.push_back (oldGranCladePath);
    mergeComponents.push_back (newBranchPath);
    const AlignPath newPath = alignPathMerge (mergeComponents);

    const GuideAlignmentEnvelope oldSibEnv = sampler.makeGuide (history.tree, nodeClosestLeaf, oldSibClosestLeaf, newPath, node, oldSibling);
    const SiblingMatrix oldSibMatrix (sampler.model, nodeSeq, oldSibSeq, parentNodeDist, parentOldSibDist, oldSibEnv, nodeEnvPos, oldSibEnvPos, node, oldSibling, parent);
    const LogProb logPostOldSiblingPath = oldSibMatrix.logPostProb (oldSiblingPath);

    const GuideAlignmentEnvelope oldBranchEnv = sampler.makeGuide (history.tree, oldGranClosestLeaf, oldParentClosestLeaf, newPath, oldGrandparent, oldParentClosestChild);
    const PosWeightMatrix oldParentSeq = oldSibMatrix.parentSeq (oldSiblingPath);
    const BranchMatrix oldBranchMatrix (sampler.model, oldGranSeq, oldParentSeq, oldGranParentDist, oldBranchEnv, oldGranEnvPos, oldParentEnvPos, oldGrandparent, parent);
    const LogProb logPostOldBranchPath = oldBranchMatrix.logPostProb (oldBranchPath);

    logForwardProposal = logFwdSibSelect + logPostNewSiblingPath + logPostNewBranchPath;
    logReverseProposal = logRevSibSelect + logPostOldSiblingPath + logPostOldBranchPath;

    LogThisAt(6,"log(Q_fwd) = " << setw(10) << logFwdSibSelect << " (#newSibling) + " << setw(10) << logPostNewSiblingPath << " (parent:node:newSibling) + " << setw(10) << logPostNewBranchPath << " (newGrandparent:parent) = " << logForwardProposal << endl);
    LogThisAt(6,"log(Q_rev) = " << setw(10) << logRevSibSelect << " (#oldSibling) + " << setw(10) << logPostOldSiblingPath << " (parent:node:oldSibling) + " << setw(10) << logPostOldBranchPath << " (oldGrandparent:parent) = " << logReverseProposal << endl);

    const vguard<FastSeq>& oldUngapped = oldAlign.ungapped;
    vguard<FastSeq> newUngapped = oldUngapped;
    if (sampler.sampleAncestralSeqs) {
      newUngapped[parent].seq = sampler.sampleSeq (newParentSeq, generator);
      logForwardProposal += sampler.logSeqPostProb (newUngapped[parent].seq, newParentSeq);
      logReverseProposal += sampler.logSeqPostProb (oldUngapped[parent].seq, oldParentSeq);
    } else
      newUngapped[parent].seq = string (alignPathResiduesInRow (newSiblingPath.at (parent)), Alignment::wildcardChar);
  
    initNewHistory (newTree, newUngapped, newPath);
  }

  // we need...
  //  newGrandparent > parent
  //  parent > newSibling
  //  parent > node
  if (parent < newSibling || parent > newGrandparent)
    newHistory = newHistory.reorder (newHistory.tree.postorderSort());

  initRatio (sampler);
}

Sampler::NodeHeightMove::NodeHeightMove (const History& history, LogProb oldLogLikelihood, const Sampler& sampler, random_engine& generator)
  : Move (NodeHeight, history, oldLogLikelihood, sampler.name)
{
  Tree newTree = history.tree;
  logForwardProposal = logReverseProposal = logJacobian = 0;

  node = Sampler::randomInternalNode (newTree, generator);
  Assert (newTree.nChildren(node) == 2, "Tree is not binary");
  leftChild = newTree.getChild (node, 0);
  rightChild = newTree.getChild (node, 1);
  parent = newTree.parentNode (node);

  LogThisAt(4,"Proposing node-height move at...\n        node #" << node << ": " << newTree.seqName(node) << "\n  left-child #" << leftChild << ": " << newTree.seqName(leftChild) << "\n right-child #" << rightChild << ": " << newTree.seqName(rightChild) << (parent >= 0 ? (string("\n      parent #") + to_string(parent) + ": " + newTree.seqName(parent)) : string()) << endl);

  const TreeBranchLength lChildDist = newTree.branchLength(leftChild);
  const TreeBranchLength rChildDist = newTree.branchLength(rightChild);
  const TreeBranchLength minChildDist = min (lChildDist, rChildDist);

  if (parent < 0) {
    const double maxLogMultiplier = log(2);
    uniform_real_distribution<double> distribution (-maxLogMultiplier, maxLogMultiplier);
    const double logMultiplier = distribution (generator);
    const double multiplier = exp (logMultiplier);
    const double newMinChildDist = minChildDist * multiplier;
    
    newTree.node[leftChild].d = lChildDist - minChildDist + newMinChildDist;
    newTree.node[rightChild].d = rChildDist - minChildDist + newMinChildDist;

    logJacobian += logMultiplier;

    LogThisAt(6,"Sampled coalescence time of #" << leftChild << " and #" << rightChild << ": " << newMinChildDist << " (previously " << minChildDist << ", multiplier " << multiplier << ")" << endl);

  } else {
    const TreeBranchLength pDist = max (0., newTree.branchLength(node));
    const TreeBranchLength pDistRange = pDist + minChildDist;

    uniform_real_distribution<TreeBranchLength> distribution (0, pDistRange);
    const TreeBranchLength pDistNew = distribution (generator);
    const TreeBranchLength cDistNew = pDistRange - pDistNew;

    newTree.node[node].d = pDistNew;
    newTree.node[leftChild].d = (lChildDist - minChildDist) + cDistNew;
    newTree.node[rightChild].d = (rChildDist - minChildDist) + cDistNew;

    LogThisAt(6,"Sampled coalescence time of #" << leftChild << " and #" << rightChild << ": " << cDistNew << " (previously " << minChildDist << ", maximum " << pDist << ")" << endl);
  }

  initNewHistory (newTree, oldHistory.gapped);
  initRatio (sampler);
}

Sampler::RescaleMove::RescaleMove (const History& history, LogProb oldLogLikelihood, const Sampler& sampler, random_engine& generator)
  : Move (Rescale, history, oldLogLikelihood, sampler.name)
{
  LogThisAt(4,"Proposing rescale move" << endl);

  const auto dist = history.tree.distanceFromRoot();
  const TreeBranchLength oldTreeHeight = *max_element (dist.begin(), dist.end());

  const double maxLogMultiplier = log(2);
  uniform_real_distribution<double> distribution (-maxLogMultiplier, maxLogMultiplier);
  const double logMultiplier = distribution (generator);
  const double multiplier = exp (logMultiplier);
  const TreeBranchLength newTreeHeight = oldTreeHeight * multiplier;

  LogThisAt(6,"Sampled tree height: " << newTreeHeight << " (previously " << oldTreeHeight << ", multiplier " << multiplier << ")" << endl);

  Tree newTree = history.tree;
  for (auto& node : newTree.node)
    node.d *= multiplier;

  logForwardProposal = logReverseProposal = 0;
  logJacobian = logMultiplier;

  initNewHistory (newTree, oldHistory.gapped);
  initRatio (sampler);
}

TreeAlignFuncs::BranchMatrixBase::BranchMatrixBase (const RateModel& model, const PosWeightMatrix& xSeq, const PosWeightMatrix& ySeq, TreeBranchLength dist, const GuideAlignmentEnvelope& env, const vguard<SeqIdx>& xEnvPos, const vguard<SeqIdx>& yEnvPos, AlignRowIndex x, AlignRowIndex y)
  : SparseDPMatrix (env, xEnvPos, yEnvPos),
    model (model),
    probModel (model, max (Tree::minBranchLength, dist)),
    logProbModel (probModel),
    xRow (x),
    yRow (y),
    xSeq (xSeq),
    ySub (Sampler::preMultiply (ySeq, logProbModel.logSubProb)),
    yEmit (Sampler::calcInsProbs (ySeq, logProbModel.logInsProb, logProbModel.logCptWeight))
{
  mm = lpTrans (ProbModel::Match, ProbModel::Match);
  mi = lpTrans (ProbModel::Match, ProbModel::Insert);
  md = lpTrans (ProbModel::Match, ProbModel::Delete);
  me = lpTrans (ProbModel::Match, ProbModel::End);

  im = lpTrans (ProbModel::Insert, ProbModel::Match);
  ii = lpTrans (ProbModel::Insert, ProbModel::Insert);
  id = lpTrans (ProbModel::Insert, ProbModel::Delete);
  ie = lpTrans (ProbModel::Insert, ProbModel::End);

  dm = lpTrans (ProbModel::Delete, ProbModel::Match);
  dd = lpTrans (ProbModel::Delete, ProbModel::Delete);
  de = lpTrans (ProbModel::Delete, ProbModel::End);

  LogThisAt(8,"Parent profile:\n" << Sampler::profileToString(xSeq)
	    << "Child profile:\n" << Sampler::profileToString(ySeq));
}

Sampler::BranchMatrix::BranchMatrix (const RateModel& model, const PosWeightMatrix& xSeq, const PosWeightMatrix& ySeq, TreeBranchLength dist, const GuideAlignmentEnvelope& env, const vguard<SeqIdx>& xEnvPos, const vguard<SeqIdx>& yEnvPos, AlignRowIndex x, AlignRowIndex y)
  : BranchMatrixBase (model, xSeq, ySeq, dist, env, xEnvPos, yEnvPos, x, y)
{
  ProgressLog (plog, 5);
  plog.initProgress ("Branch alignment matrix (%u*%u)", xSize, ySize);

  lpStart() = 0;
  for (SeqIdx xpos = 0; xpos < xSize; ++xpos) {

    plog.logProgress (xpos / (double) (xSize - 1), "row %d/%d", xpos + 1, xSize);

    for (SeqIdx ypos = 0; ypos < ySize; ++ypos)
      if (inEnvelope (xpos, ypos)) {
	XYCell& dest = xyCell (xpos, ypos);

	if (xpos > 0 && inEnvelope (xpos - 1, ypos)) {
	  const XYCell& dSrc = xyCell (xpos - 1, ypos);

	  dest(ProbModel::Delete) = log_sum_exp (dSrc(ProbModel::Match) + md,
						 dSrc(ProbModel::Insert) + id,
						 dSrc(ProbModel::Delete) + dd);
	}

      	if (ypos > 0 && inEnvelope (xpos, ypos - 1)) {
	  const XYCell& iSrc = xyCell (xpos, ypos - 1);
	  const LogProb yEmitScore = yEmit[ypos - 1];

	  dest(ProbModel::Insert) = yEmitScore + log_sum_exp (iSrc(ProbModel::Match) + mi,
							      iSrc(ProbModel::Insert) + ii);
	}

	if (xpos > 0 && ypos > 0 && inEnvelope (xpos - 1, ypos - 1)) {
	  const XYCell& mSrc = xyCell (xpos - 1, ypos - 1);
	  const LogProb xyEmitScore = logMatch (xpos, ypos);
	  
	  dest(ProbModel::Match) = xyEmitScore + log_sum_exp (mSrc(ProbModel::Match) + mm,
							      mSrc(ProbModel::Insert) + im,
							      mSrc(ProbModel::Delete) + dm);
	}
      }
  }

  const XYCell& endCell = xyCell (xSize - 1, ySize - 1);

  lpEnd = log_sum_exp (endCell(ProbModel::Match) + me,
		       endCell(ProbModel::Insert) + ie,
		       endCell(ProbModel::Delete) + de);

  if (LoggingThisAt(9))
    writeToLog(9);

  LogThisAt(6,"Forward log-likelihood is " << lpEnd << endl);
}

AlignPath Sampler::BranchMatrix::sample (random_engine& generator) const {
  CellCoords coords (xSize - 1, ySize - 1, ProbModel::End);
  AlignRowPath xPath, yPath;
  while (coords.xpos > 0 || coords.ypos > 0) {
    bool x, y;
    getColumn (coords, x, y);
    if (x || y) {
      xPath.push_back (x);
      yPath.push_back (y);
    }
    CellCoords src = coords;
    if (x) --src.xpos;
    if (y) --src.ypos;
    const XYCell& srcCell = xyCell (src.xpos, src.ypos);
    const LogProb e = lpEmit (coords);
    map<CellCoords,LogProb> srcLogProb;
    for (src.state = 0; src.state < (unsigned int) ProbModel::End; ++src.state)
      srcLogProb[src] = srcCell(src.state) + lpTrans ((State) src.state, (State) coords.state) + e;

    const double lpTot = log_sum_exp (extract_values (srcLogProb));
    Assert (SAMPLER_NEAR_EQ (lpTot, cell(coords)), "Traceback total (%g) doesn't match stored value (%g) at cell %s", lpTot, cell(coords), coords.toString().c_str());
    Assert (lpTot > -numeric_limits<double>::infinity(), "Traceback state has zero probability at cell %s", coords.toString().c_str());

    coords = random_key_log (srcLogProb, generator);
  }
  AlignPath path;
  path[xRow] = AlignRowPath (xPath.rbegin(), xPath.rend());
  path[yRow] = AlignRowPath (yPath.rbegin(), yPath.rend());

  LogThisAt(9,"Sampled parent-child alignment:" << endl << alignPathString(path));

  return path;
}

LogProb Sampler::BranchMatrixBase::logPathProb (const AlignPath& path) const {
  const AlignColIndex cols = alignPathColumns (path);
  LogProb lp = 0, lpPathTrans = 0, lpPathEmit = 0;
  CellCoords c (0, 0, ProbModel::Start);
  vguard<LogProb> colEmit (cols, 0);
  for (AlignColIndex col = 0; col < cols; ++col) {
    const bool dx = path.at(xRow)[col];
    const bool dy = path.at(yRow)[col];
    if (dx)
      ++c.xpos;
    if (dy)
      ++c.ypos;
    const ProbModel::State prevState = (ProbModel::State) c.state;
    c.state = ProbModel::getState (dx, dy);
    if (!inEnvelope (c.xpos, c.ypos))
      return -numeric_limits<double>::infinity();
    const LogProb lpt = lpTrans (prevState, (ProbModel::State) c.state);
    const LogProb lpe = lpEmit (c);
    lp += lpt + lpe;
    Assert (lp <= cell(c) * (1 - SAMPLER_EPSILON), "Positive posterior probability");
    lp = min (lp, cell(c));  // mitigate precision errors
    lpPathTrans += lpt;
    lpPathEmit += lpe;
    colEmit[col] = lpe;
  }
  const LogProb lpte = lpTrans ((ProbModel::State) c.state, ProbModel::End);
  lp += lpte;
  lpPathTrans += lpte;
  LogThisAt(8,"Path:\n" << alignPathString(path)
	    << "Log-likelihood = " << lp << " = " << lpPathTrans << " (trans) + " << lpPathEmit << " (emit)\n"
	    << "Column emissions: (" << to_string_join(colEmit) << ")" << endl);
  return lp;
}

LogProb Sampler::BranchMatrix::logPostProb (const AlignPath& path) const {
  const LogProb lp = logPathProb (path);
  Assert (lp <= lpEnd * (1 - SAMPLER_EPSILON), "Positive posterior probability");
  return min (lp, lpEnd) - lpEnd;
}

LogProb TreeAlignFuncs::BranchMatrixBase::lpTrans (State src, State dest) const {
  return log (probModel.transProb (src, dest));
}

LogProb TreeAlignFuncs::BranchMatrixBase::lpEmit (const CellCoords& coords) const {
  switch ((State) coords.state) {
  case ProbModel::Match: return coords.xpos > 0 && coords.ypos > 0 ? logMatch (coords.xpos, coords.ypos) : -numeric_limits<double>::infinity();
  case ProbModel::Insert: return coords.ypos > 0 ? yEmit[coords.ypos - 1] : -numeric_limits<double>::infinity();
  default: break;
  }
  return 0;
}

void TreeAlignFuncs::BranchMatrixBase::getColumn (const CellCoords& coords, bool& x, bool& y) {
  x = y = false;
  switch ((State) coords.state) {
  case ProbModel::Match: if (coords.xpos > 0 && coords.ypos > 0) x = y = true; break;
  case ProbModel::Insert: y = true; break;
  case ProbModel::Delete: x = true; break;
  default: break;
  }
}

Sampler::SiblingMatrix::SiblingMatrix (const RateModel& model, const PosWeightMatrix& lSeq, const PosWeightMatrix& rSeq, TreeBranchLength plDist, TreeBranchLength prDist, const GuideAlignmentEnvelope& env, const vguard<SeqIdx>& lEnvPos, const vguard<SeqIdx>& rEnvPos, AlignRowIndex l, AlignRowIndex r, AlignRowIndex p)
  : SparseDPMatrix (env, lEnvPos, rEnvPos),
    model (model),
    lProbModel (model, max (Tree::minBranchLength, plDist)),
    rProbModel (model, max (Tree::minBranchLength, prDist)),
    lLogProbModel (lProbModel),
    rLogProbModel (rProbModel),
    logRoot (log_vector_gsl_vector (model.insProb)),
    lRow (l),
    rRow (r),
    pRow (p),
    lSub (Sampler::preMultiply (lSeq, lLogProbModel.logSubProb)),
    rSub (Sampler::preMultiply (rSeq, rLogProbModel.logSubProb)),
    lEmit (Sampler::calcInsProbs (lSeq, lLogProbModel.logInsProb, lLogProbModel.logCptWeight)),
    rEmit (Sampler::calcInsProbs (rSeq, rLogProbModel.logInsProb, rLogProbModel.logCptWeight))
{
  for (int cpt = 0; cpt < model.components(); ++cpt)
    for (AlphTok tok = 0; tok < model.alphabetSize(); ++tok)
      ((vguard<vguard<LogProb> >&)logRoot)[cpt][tok] += log (model.cptWeight[cpt]);

  imm_www = lpTransElimSelfLoopIDD (IMM, WWW);
  imm_imi = lpTransElimSelfLoopIDD (IMM, IMI);
  imm_iiw = lpTransElimSelfLoopIDD (IMM, IIW);

  imd_wwx = lpTransElimSelfLoopIDD (IMD, WWX);
  imd_iix = lpTransElimSelfLoopIDD (IMD, IIX);

  idm_wxw = lpTransElimSelfLoopIDD (IDM, WXW);
  idm_idi = lpTransElimSelfLoopIDD (IDM, IDI);

  idd_imm = lpTransElimSelfLoopIDD (IDD, IMM);
  idd_imd = lpTransElimSelfLoopIDD (IDD, IMD);
  idd_idm = lpTransElimSelfLoopIDD (IDD, IDM);
  idd_eee = lpTransElimSelfLoopIDD (IDD, EEE);

  www_imm = lpTransElimSelfLoopIDD (WWW, IMM);
  www_imd = lpTransElimSelfLoopIDD (WWW, IMD);
  www_idm = lpTransElimSelfLoopIDD (WWW, IDM);
  www_idd = lpTransElimSelfLoopIDD (WWW, IDD);
  www_eee = lpTransElimSelfLoopIDD (WWW, EEE);

  wwx_imm = lpTransElimSelfLoopIDD (WWX, IMM);
  wwx_imd = lpTransElimSelfLoopIDD (WWX, IMD);
  wwx_idm = lpTransElimSelfLoopIDD (WWX, IDM);
  wwx_idd = lpTransElimSelfLoopIDD (WWX, IDD);
  wwx_eee = lpTransElimSelfLoopIDD (WWX, EEE);

  wxw_imm = lpTransElimSelfLoopIDD (WXW, IMM);
  wxw_imd = lpTransElimSelfLoopIDD (WXW, IMD);
  wxw_idm = lpTransElimSelfLoopIDD (WXW, IDM);
  wxw_idd = lpTransElimSelfLoopIDD (WXW, IDD);
  wxw_eee = lpTransElimSelfLoopIDD (WXW, EEE);

  imi_www = lpTransElimSelfLoopIDD (IMI, WWW);
  imi_imi = lpTransElimSelfLoopIDD (IMI, IMI);
  imi_iiw = lpTransElimSelfLoopIDD (IMI, IIW);

  iiw_www = lpTransElimSelfLoopIDD (IIW, WWW);
  iiw_iiw = lpTransElimSelfLoopIDD (IIW, IIW);

  idi_wxw = lpTransElimSelfLoopIDD (IDI, WXW);
  idi_idi = lpTransElimSelfLoopIDD (IDI, IDI);

  iix_wwx = lpTransElimSelfLoopIDD (IIX, WWX);
  iix_iix = lpTransElimSelfLoopIDD (IIX, IIX);

  LogThisAt(8,"Left-node profile:\n" << Sampler::profileToString(lSeq)
	    << "Right-node profile:\n" << Sampler::profileToString(rSeq));

  ProgressLog (plog, 5);
  plog.initProgress ("Parent proposal matrix (%u*%u)", xSize, ySize);

  lpStart() = 0;
  cell(0,0,WWW) = imm_www;
  for (SeqIdx xpos = 0; xpos < xSize; ++xpos) {

    plog.logProgress (xpos / (double) (xSize - 1), "row %d/%d", xpos + 1, xSize);

    for (SeqIdx ypos = 0; ypos < ySize; ++ypos)
      if (inEnvelope (xpos, ypos)) {
	XYCell& dest = xyCell (xpos, ypos);

	if (xpos > 0 && inEnvelope (xpos - 1, ypos)) {
	  const XYCell& lSrc = xyCell (xpos - 1, ypos);
	  const LogProb lEmitScore = lEmit[xpos - 1];

	  dest(IIW) = lEmitScore + log_sum_exp (lSrc(IMM) + imm_iiw,
						lSrc(IMI) + imi_iiw,
						lSrc(IIW) + iiw_iiw);

	  dest(IIX) = lEmitScore + log_sum_exp (lSrc(IMD) + imd_iix,
						lSrc(IIX) + iix_iix);

	  dest(IMD) = lEmitScore + log_sum_exp (lSrc(WWW) + www_imd,
						lSrc(WWX) + wwx_imd,
						lSrc(WXW) + wxw_imd,
						lSrc(IDD) + idd_imd);

	  dest(WWW) = dest(IIW) + iiw_www;

	  dest(WWX) = log_sum_exp (dest(IIX) + iix_wwx,
				   dest(IMD) + imd_wwx);
	}

      	if (ypos > 0 && inEnvelope (xpos, ypos - 1)) {
	  const XYCell& rSrc = xyCell (xpos, ypos - 1);
	  const LogProb rEmitScore = rEmit[ypos - 1];

	  dest(IMI) = rEmitScore + log_sum_exp (rSrc(IMM) + imm_imi,
						rSrc(IMI) + imi_imi);

	  dest(IDI) = rEmitScore + log_sum_exp (rSrc(IDM) + idm_idi,
						rSrc(IDI) + idi_idi);

	  dest(IDM) = rEmitScore + log_sum_exp (rSrc(WWW) + www_idm,
						rSrc(WWX) + wwx_idm,
						rSrc(WXW) + wxw_idm,
						rSrc(IDD) + idd_idm);

	  dest(WWW) = log_sum_exp (dest(WWW),
				   dest(IMI) + imi_www);

	  dest(WXW) = log_sum_exp (dest(IDI) + idi_wxw,
				   dest(IDM) + idm_wxw);
	}

	if (xpos > 0 && ypos > 0 && inEnvelope (xpos - 1, ypos - 1)) {
	  const XYCell& lrSrc = xyCell (xpos - 1, ypos - 1);
	  const LogProb lrEmitScore = logMatch (xpos, ypos);
	  
	  dest(IMM) = lrEmitScore + log_sum_exp (lrSrc(WWW) + www_imm,
						 lrSrc(WWX) + wwx_imm,
						 lrSrc(WXW) + wxw_imm,
						 lrSrc(IDD) + idd_imm);

	  dest(WWW) = log_sum_exp (dest(WWW),
				   dest(IMM) + imm_www);
	}

	dest(IDD) = log_sum_exp (dest(WWW) + www_idd,
				 dest(WWX) + wwx_idd,
				 dest(WXW) + wxw_idd);
      }
  }

  const XYCell& endCell = xyCell (xSize - 1, ySize - 1);

  lpEnd = log_sum_exp (endCell(IDD) + idd_eee,
		       endCell(WWW) + www_eee,
		       endCell(WWX) + wwx_eee,
		       endCell(WXW) + wxw_eee);

  if (LoggingThisAt(9))
    writeToLog(9);

  LogThisAt(6,"Forward log-likelihood is " << lpEnd << endl);
}

AlignPath Sampler::SiblingMatrix::sample (random_engine& generator) const {
  CellCoords coords (xSize - 1, ySize - 1, EEE);
  AlignRowPath lPath, rPath, pPath;
  while (coords.xpos > 0 || coords.ypos > 0) {
    bool l, r, p;
    getColumn (coords, l, r, p);
    if (l || r || p) {
      lPath.push_back (l);
      rPath.push_back (r);
      pPath.push_back (p);
    }
    if ((State) coords.state == IDD) {  // explicitly add IDD self-loops outside main traceback step
      geometric_distribution<int> distribution (iddSelfLoopProb());
      int iddSelfLoops = distribution (generator);
      while (iddSelfLoops-- > 0) {
	lPath.push_back (l);
	rPath.push_back (r);
	pPath.push_back (p);
      }
    }
    CellCoords src = coords;
    if (l) --src.xpos;
    if (r) --src.ypos;
    const XYCell& srcCell = xyCell (src.xpos, src.ypos);
    const LogProb e = lpEmit (coords);
    map<CellCoords,LogProb> srcLogProb;
    for (src.state = 0; src.state < (unsigned int) EEE; ++src.state)
      srcLogProb[src] = srcCell(src.state) + lpTransElimSelfLoopIDD ((State) src.state, (State) coords.state) + e;

    const double lpTot = log_sum_exp (extract_values (srcLogProb));
    Assert (SAMPLER_NEAR_EQ (lpTot, cell(coords)), "Traceback total (%g) doesn't match stored value (%g) at cell %s", lpTot, cell(coords), coords.toString().c_str());
    Assert (lpTot > -numeric_limits<double>::infinity(), "Traceback state has zero probability at cell %s", coords.toString().c_str());

    coords = random_key_log (srcLogProb, generator);
  }
  AlignPath path;
  path[lRow] = AlignRowPath (lPath.rbegin(), lPath.rend());
  path[rRow] = AlignRowPath (rPath.rbegin(), rPath.rend());
  path[pRow] = AlignRowPath (pPath.rbegin(), pPath.rend());

  LogThisAt(9,"Sampled sibling-parent alignment:" << endl << alignPathString(path));

  return path;
}

LogProb Sampler::SiblingMatrix::logPostProb (const AlignPath& lrpPath) const {
  const AlignColIndex cols = alignPathColumns (lrpPath);
  LogProb lp = 0;
  CellCoords c (0, 0, SSS);
  for (AlignColIndex col = 0; col < cols; ++col) {
    const bool dl = lrpPath.at(lRow)[col];
    const bool dr = lrpPath.at(rRow)[col];
    const bool dp = lrpPath.at(pRow)[col];
    if (dl)
      ++c.xpos;
    if (dr)
      ++c.ypos;
    const State prevState = (State) c.state;
    c.state = getState (prevState, dl, dr, dp);
    if (!inEnvelope (c.xpos, c.ypos))
      return -numeric_limits<double>::infinity();
    lp += lpTransElimWait (prevState, (State) c.state) + lpEmit (c);
    Assert (lp <= cell(c) * (1 - SAMPLER_EPSILON), "Positive posterior probability");
    lp = min (lp, cell(c));  // mitigate precision errors
  }
  lp += lpTransElimWait ((State) c.state, EEE);
  Assert (lp <= lpEnd * (1 - SAMPLER_EPSILON), "Positive posterior probability");
  lp = min (lp, lpEnd);
  return lp - lpEnd;
}

LogProb Sampler::SiblingMatrix::lpEmit (const CellCoords& coords) const {
  switch ((State) coords.state) {
  case IMM: return coords.xpos > 0 && coords.ypos > 0 ? logMatch (coords.xpos, coords.ypos) : -numeric_limits<double>::infinity();
  case IDM: case IMI: case IDI: return coords.ypos > 0 ? rEmit[coords.ypos - 1] : -numeric_limits<double>::infinity();
  case IMD: case IIW: case IIX: return coords.xpos > 0 ? lEmit[coords.xpos - 1] : -numeric_limits<double>::infinity();
  default: break;
  }
  return 0;
}

Sampler::SiblingMatrix::State Sampler::SiblingMatrix::getState (State src, bool leftUngapped, bool rightUngapped, bool parentUngapped) {
  if (parentUngapped)
    return leftUngapped ? (rightUngapped ? IMM : IMD) : (rightUngapped ? IDM : IDD);
  if (leftUngapped)
    return (src == IMD || src == IIX) ? IIX : IIW;
  if (rightUngapped)
    return (src == IDM || src == IDI) ? IDI : IMI;
  if (src == IDM || src == IDD || src == IDI)
    return WXW;
  if (src == IMD || src == IIX)
    return WWX;
  return WWW;
}

void Sampler::SiblingMatrix::getColumn (const CellCoords& coords, bool& l, bool& r, bool& p) {
  p = l = r = false;
  switch ((State) coords.state) {
  case IMM: if (coords.xpos > 0 && coords.ypos > 0) p = l = r = true; break;
  case IMD: p = l = true; break;
  case IDM: p = r = true; break;
  case IDD: p = true; break;
  case IIW: case IIX: if (coords.xpos > 0) l = true; break;
  case IMI: case IDI: if (coords.ypos > 0) r = true; break;
  default: break;
  }
}

LogProb Sampler::SiblingMatrix::lpTransElimSelfLoopIDD (State src, State dest) const {
  return src == IDD
    ? (dest == IDD
       ? -numeric_limits<double>::infinity()
       : lpTrans(src,dest) + iddExit())
    : lpTrans(src,dest);
}

LogProb Sampler::SiblingMatrix::lpTrans (State src, State dest) const {
  switch (src) {
  case IMM:
    switch (dest) {
    case WWW: return lNoIns() + rNoIns();
    case IMI: return rIns();
    case IIW: return lIns() + rNoIns();
    default: break;
    }
    break;

  case IMD:
    switch (dest) {
    case WWX: return lNoIns();
    case IIX: return lIns();
    default: break;
    }
    break;

  case IDM:
    switch (dest) {
    case WXW: return rNoIns();
    case IDI: return rIns();
    default: break;
    }
    break;

  case IDD:
    switch (dest) {
    case IDD: return iddStay();
    case IMM: return rootExt() + lNoDelExt() + rNoDelExt();
    case IMD: return rootExt() + lNoDelExt() + rDelExt();
    case IDM: return rootExt() + lDelExt() + rNoDelExt();
    case EEE: return rootNoExt() + lNoDelExt() + rNoDelExt();
    default: break;
    }
    break;

  case WWW:
    switch (dest) {
    case IMM: return rootExt() + lNoDel() + rNoDel();
    case IMD: return rootExt() + lNoDel() + rDel();
    case IDM: return rootExt() + lDel() + rNoDel();
    case IDD: return rootExt() + lDel() + rDel();
    case EEE: return 0;
    default: break;
    }
    break;

  case WWX:
    switch (dest) {
    case IMM: return rootExt() + lNoDel() + rNoDelExt();
    case IMD: return rootExt() + lNoDel() + rDelExt();
    case IDM: return rootExt() + lDel() + rNoDelExt();
    case IDD: return rootExt() + lDel() + rDelExt();
    case EEE: return rNoDelExt();
    default: break;
    }
    break;

  case WXW:
    switch (dest) {
    case IMM: return rootExt() + lNoDelExt() + rNoDel();
    case IMD: return rootExt() + lNoDelExt() + rDel();
    case IDM: return rootExt() + lDelExt() + rNoDel();
    case IDD: return rootExt() + lDelExt() + rDel();
    case EEE: return lNoDelExt();
    default: break;
    }
    break;

  case IMI:
    switch (dest) {
    case WWW: return lNoIns() + rNoInsExt();
    case IMI: return rInsExt();
    case IIW: return lIns() + rNoInsExt();
    default: break;
    }
    break;

  case IIW:
    switch (dest) {
    case WWW: return lNoInsExt();
    case IIW: return lInsExt();
    default: break;
    }
    break;

  case IDI:
    switch (dest) {
    case WXW: return rNoInsExt();
    case IDI: return rInsExt();
    default: break;
    }
    break;

  case IIX:
    switch (dest) {
    case WWX: return lNoInsExt();
    case IIX: return lInsExt();
    default: break;
    }
    break;

  default:
    break;
  }
  return -numeric_limits<double>::infinity();
}

LogProb Sampler::SiblingMatrix::lpTransElimWait (State src, State dest) const {
  return log_sum_exp (lpTrans (src, dest),
		      lpTrans (src, WWW) + lpTrans (WWW, dest),
		      lpTrans (src, WWX) + lpTrans (WWX, dest),
		      lpTrans (src, WXW) + lpTrans (WXW, dest));
}

Sampler::PosWeightMatrix Sampler::SiblingMatrix::parentSeq (const AlignPath& lrpPath) const {
  PosWeightMatrix pwm;
  pwm.reserve (alignPathResiduesInRow (lrpPath.at (pRow)));
  const AlignColIndex cols = alignPathColumns (lrpPath);
  SeqIdx lPos = 0, rPos = 0;
  for (AlignColIndex col = 0; col < cols; ++col)
    if (lrpPath.at(pRow)[col]) {
      vguard<vguard<LogProb> > prof (model.components(), vguard<LogProb> (model.alphabetSize(), 0));
      if (lrpPath.at(lRow)[col]) {
	for (int cpt = 0; cpt < model.components(); ++cpt)
	  for (AlphTok i = 0; i < model.alphabetSize(); ++i)
	    prof[cpt][i] += lSub[lPos][cpt][i];
	++lPos;
      }
      if (lrpPath.at(rRow)[col]) {
	for (int cpt = 0; cpt < model.components(); ++cpt)
	  for (AlphTok i = 0; i < model.alphabetSize(); ++i)
	    prof[cpt][i] += rSub[rPos][cpt][i];
	++rPos;
      }
      LogProb norm = -numeric_limits<double>::infinity();
      for (int cpt = 0; cpt < model.components(); ++cpt)
	for (AlphTok i = 0; i < model.alphabetSize(); ++i)
	  log_accum_exp (norm, prof[cpt][i]);
      for (int cpt = 0; cpt < model.components(); ++cpt)
	for (AlphTok i = 0; i < model.alphabetSize(); ++i)
	  prof[cpt][i] -= norm;
      pwm.push_back (prof);
    }
  return pwm;
}

Sampler::Sampler (const RateModel& model, const SimpleTreePrior& treePrior, const vguard<FastSeq>& gappedGuide)
  : model (model),
    treePrior (treePrior),
    moveRate (Move::TotalMoveTypes, 1.),
    movesProposed (Move::TotalMoveTypes, 0),
    movesAccepted (Move::TotalMoveTypes, 0),
    moveNanosecs (Move::TotalMoveTypes, 0.),
    useFixedGuide (false),
    sampleAncestralSeqs (false),
    guide (gappedGuide),
    maxDistanceFromGuide (DefaultMaxDistanceFromGuide)
{
  for (AlignRowIndex r = 0; r < gappedGuide.size(); ++r) {
    const string& name = gappedGuide[r].name;
    Assert (guideRowByName.count(name) == 0, "Duplicate name %s in guide alignment", name.c_str());
    guideRowByName[name] = r;
  }
}

Sampler::Move Sampler::proposeMove (const History& oldHistory, LogProb oldLogLikelihood, random_engine& generator) const {
  const Move::Type type = (Move::Type) random_index (moveRate, generator);
  switch (type) {
  case Move::BranchAlign: return BranchAlignMove (oldHistory, oldLogLikelihood, *this, generator);
  case Move::NodeAlign: return NodeAlignMove (oldHistory, oldLogLikelihood, *this, generator);
  case Move::PruneAndRegraft: return PruneAndRegraftMove (oldHistory, oldLogLikelihood, *this, generator);
  case Move::NodeHeight: return NodeHeightMove (oldHistory, oldLogLikelihood, *this, generator);
  case Move::Rescale: return RescaleMove (oldHistory, oldLogLikelihood, *this, generator);
  default: break;
  }
  Abort ("Unknown move type");
  return Move();
}

void Sampler::addLogger (Logger& logger) {
  loggers.push_back (&logger);
}

void Sampler::initialize (const History& initialHistory, const string& samplerName) {
  name = samplerName;
  currentHistory = initialHistory;
  currentHistory.assertNamesMatch();

  isUltrametric = currentHistory.tree.isUltrametric();
  if (isUltrametric)
    LogThisAt(3,"Initial tree is ultrametric" << endl);
  else
    LogThisAt(1,"WARNING: initial tree is not ultrametric" << endl);
  
  bestHistory = currentHistory;
  currentLogLikelihood = bestLogLikelihood = logLikelihood (currentHistory, "initial");

  // set move rates more-or-less arbitrarily
  moveRate[Move::BranchAlign] = initialHistory.tree.hasChildren() ? 1 : 0;
  moveRate[Move::NodeAlign] = 1;
  moveRate[Move::PruneAndRegraft] = initialHistory.tree.hasGrandchildren() ? 1 : 0;
  moveRate[Move::NodeHeight] = 2;
  moveRate[Move::Rescale] = 2;
}

void Sampler::fixTree() {
  moveRate[Move::PruneAndRegraft] = 0;
  moveRate[Move::NodeHeight] = 0;
  moveRate[Move::Rescale] = 0;
}

void Sampler::fixAlignment() {
  moveRate[Move::BranchAlign] = 0;
  moveRate[Move::NodeAlign] = 0;
}

void Sampler::sample (random_engine& generator) {
    // propose
    const std::chrono::system_clock::time_point before = std::chrono::system_clock::now();
    const Move move = proposeMove (currentHistory, currentLogLikelihood, generator);
    const std::chrono::system_clock::time_point after = std::chrono::system_clock::now();
    moveNanosecs[move.type] += std::chrono::duration_cast<std::chrono::nanoseconds> (after - before).count();
    ++movesProposed[move.type];

    // do some consistency checks
    move.newHistory.assertNamesMatch();
    move.newHistory.tree.assertPostorderSorted();
    if (isUltrametric && !move.newHistory.tree.isUltrametric())
      Warn ("Move generated a non-ultrametric tree");
    
    // accept/reject
    if (move.accept (generator)) {
      currentHistory = move.newHistory;
      currentLogLikelihood = move.newLogLikelihood;
      ++movesAccepted[move.type];
    }

    // log
    for (auto& logger : loggers)
      logger->logHistory (currentHistory);

    // keep track of best history
    if (move.newLogLikelihood > bestLogLikelihood) {
      bestHistory = move.newHistory;
      bestLogLikelihood = move.newLogLikelihood;
      LogThisAt(2,"New best log-likelihood: " << bestLogLikelihood << " (" << name << ")" << endl);
    }
}

void Sampler::run (vguard<Sampler>& samplers, random_engine& generator, unsigned int nSamples) {
  ProgressLog (plog, 2);
  plog.initProgress ("MCMC sampling run");

  vguard<double> nodes;
  for (const auto& sampler: samplers)
    nodes.push_back (sampler.currentHistory.tree.nodes());
  
  for (unsigned int n = 0; n < nSamples; ++n) {
    // print progress
    plog.logProgress (n / (double) (nSamples - 1), "step %u/%u", n + 1, nSamples);

    // select a sampler, weighted by # of nodes
    const size_t nSampler = random_index (nodes, generator);
    LogThisAt(4,"Sampling dataset #" << nSampler+1 << ": " << samplers[nSampler].name << endl);
    
    // sample
    samplers[nSampler].sample (generator);
  }

  // log stats
  for (size_t nSampler = 0; nSampler < samplers.size(); ++nSampler)
    LogThisAt(1,"Dataset #" << nSampler+1 << " (" << samplers[nSampler].name << "):\n" << samplers[nSampler].moveStats());
}

string Sampler::moveStats() const {
  ostringstream out;
  for (int t = 0; t < (int) Move::TotalMoveTypes; ++t)
    out << setw(Move::typeNameWidth()) << Move::typeName ((Move::Type) t) << ": "
	<< setw(5) << movesProposed[t] << " moves, "
	<< setw(5) << movesAccepted[t] << " accepted, "
	<< setw(12) << (moveNanosecs[t] / 1e9) << " seconds, "
	<< setw(12) << ((double) movesAccepted[t] / (moveNanosecs[t] / 1e9)) << " accepted/sec"
	<< endl;
  return out.str();
}

string Sampler::sampleSeq (const PosWeightMatrix& profile, random_engine& generator) const {
  string seq (profile.size(), Alignment::wildcardChar);
  for (SeqIdx pos = 0; pos < profile.size(); ++pos) {
    const LogProb norm = log_sum_exp (profile[pos]);
    vguard<double> p (model.alphabetSize(), 0.);
    for (int cpt = 0; cpt < model.components(); ++cpt)
      for (AlphTok tok = 0; tok < model.alphabetSize(); ++tok)
	p[tok] += exp (profile[pos][cpt][tok] - norm);
    discrete_distribution<AlphTok> tokDistrib (p.begin(), p.end());
    seq[pos] = model.alphabet[tokDistrib(generator)];
  }
  return seq;
}

LogProb Sampler::logSeqPostProb (const string& seq, const PosWeightMatrix& profile) const {
  Assert (seq.size() == profile.size(), "Sequence length (%d) does not match profile (%d)", seq.size(), profile.size());
  LogProb lp = 0;
  for (SeqIdx pos = 0; pos < profile.size(); ++pos) {
    const char c = seq[pos];
    if (!Alignment::isWildcard(c)) {
      const UnvalidatedAlphTok tok = model.tokenize (c);
      if (tok < 0)
	return -numeric_limits<double>::infinity();
      const LogProb norm = log_sum_exp (profile[pos]);
      double lpTok = -numeric_limits<double>::infinity();
      for (int cpt = 0; cpt < model.components(); ++cpt)
	log_accum_exp (lpTok, profile[pos][cpt][tok] - norm);
      lp += lpTok;
    }
  }
  return lp;
}

string Sampler::profileToString (const PosWeightMatrix& profile) {
  ostringstream out;
  for (const auto& wm: profile) {
    for (int cpt = 0; cpt < (int) wm.size(); ++cpt)
      for (AlphTok i = 0; i < wm[cpt].size(); ++i)
	out << setw(10) << wm[cpt][i];
    out << endl;
  }
  return out.str();
}
