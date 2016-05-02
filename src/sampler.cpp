#include "sampler.h"
#include "util.h"

TreeNodeIndex Sampler::randomInternalNode (const Tree& tree, random_engine& generator) {
  vguard<TreeNodeIndex> intNodes;
  intNodes.reserve (tree.nodes() / 2);
  for (TreeNodeIndex n = 0; n < tree.nodes(); ++n)
    if (!tree.isLeaf(n))
      intNodes.push_back (n);
  return random_element (intNodes, generator);
}

vguard<SeqIdx> Sampler::guideSeqPos (const AlignPath& path, AlignRowIndex row, AlignRowIndex guideRow) {
  vguard<SeqIdx> guidePos;
  const auto cols = alignPathColumns (path);
  guidePos.reserve (cols);
  const AlignRowPath& rowPath = path.at(row);
  const AlignRowPath& guideRowPath = path.at(guideRow);
  SeqIdx pos = 0;
  for (AlignColIndex col = 0; col < cols; ++col) {
    if (guideRowPath[col])
      ++pos;
    if (rowPath[col])
      guidePos.push_back (pos > 0 ? (pos - 1) : 0);
  }
  guidePos.push_back (pos > 0 ? (pos - 1) : 0);
  return guidePos;
}

TokSeq Sampler::removeGapsAndTokenize (const FastSeq& gapped) const {
  TokSeq tok;
  tok.reserve (gapped.length());
  for (auto c : gapped.seq)
    if (!Alignment::isGap (c))
      tok.push_back (tokenize (c, model.alphabet));
  return tok;
}

AlignPath Sampler::cladePath (const AlignPath& path, const Tree& tree, TreeNodeIndex cladeRoot, TreeNodeIndex cladeRootParent) {
  AlignPath cp;
  for (auto n : tree.rerootedPreorderSort(cladeRoot,cladeRootParent))
    cp[n] = path.at(n);
  return cp;
}

AlignPath Sampler::treePathToSiblingPath (const AlignPath& path, const Tree& tree, TreeNodeIndex parent) {
  AlignPath p;
  Assert (tree.nChildren(parent) == 2, "Non-binary tree");
  p[SiblingMatrix::pRow()] = path.at(parent);
  p[SiblingMatrix::lRow()] = path.at(tree.getChild(parent,0));
  p[SiblingMatrix::rRow()] = path.at(tree.getChild(parent,1));
  return p;
}

AlignPath Sampler::siblingPathToTreePath (const AlignPath& path, const Tree& tree, TreeNodeIndex parent) {
  AlignPath p;
  Assert (tree.nChildren(parent) == 2, "Non-binary tree");
  p[parent] = path.at(SiblingMatrix::pRow());
  p[tree.getChild(parent,0)] = path.at(SiblingMatrix::lRow());
  p[tree.getChild(parent,1)] = path.at(SiblingMatrix::rRow());
  return p;
}

AlignPath Sampler::treePathToBranchPath (const AlignPath& path, const Tree& tree, TreeNodeIndex node) {
  AlignPath p;
  p[BranchMatrix::xRow()] = path.at(tree.parentNode(node));
  p[BranchMatrix::yRow()] = path.at(node);
  return p;
}

AlignPath Sampler::branchPathToTreePath (const AlignPath& path, const Tree& tree, TreeNodeIndex node) {
  AlignPath p;
  p[tree.parentNode(node)] = path.at(BranchMatrix::xRow());
  p[node] = path.at(BranchMatrix::yRow());
  return p;
}

Sampler::Move::Move (Type type, const History& history)
  : type (type),
    oldHistory (history)
{ }

Sampler::SampleNodeMove::SampleNodeMove (const History& history, Sampler& sampler, random_engine& generator)
  : Move (SampleNode, history)
{
  node = Sampler::randomInternalNode (history.tree, generator);
  parent = history.tree.parentNode (node);
  Assert (history.tree.nChildren(node) == 2, "Non-binary tree");
  leftChild = history.tree.getChild (node, 0);
  rightChild = history.tree.getChild (node, 1);

  const TreeBranchLength lDist = history.tree.branchLength(node,leftChild);
  const TreeBranchLength rDist = history.tree.branchLength(node,rightChild);
  
  const TreeNodeIndex leftChildClosestLeaf = history.tree.closestLeaf (leftChild, node);
  const TreeNodeIndex rightChildClosestLeaf = history.tree.closestLeaf (rightChild, node);

  const Alignment oldAlign (history.gapped);
  const vguard<SeqIdx> leftChildEnvPos = guideSeqPos (oldAlign.path, leftChild, leftChildClosestLeaf);
  const vguard<SeqIdx> rightChildEnvPos = guideSeqPos (oldAlign.path, rightChild, rightChildClosestLeaf);

  const GuideAlignmentEnvelope siblingEnv (sampler.guide.path, leftChildClosestLeaf, rightChildClosestLeaf, sampler.maxDistanceFromGuide);

  const TokSeq lSeq = sampler.removeGapsAndTokenize (history.gapped[leftChild]);
  const TokSeq rSeq = sampler.removeGapsAndTokenize (history.gapped[rightChild]);

  const AlignPath lCladePath = Sampler::cladePath (oldAlign.path, history.tree, leftChild, node);
  const AlignPath rCladePath = Sampler::cladePath (oldAlign.path, history.tree, rightChild, node);

  const SiblingMatrix sibMatrix (sampler.model, lSeq, rSeq, lDist, rDist, siblingEnv, leftChildEnvPos, rightChildEnvPos);
  const AlignPath newSiblingPath = sibMatrix.sampleAlign (generator);
  const LogProb logPostNewSiblingPath = sibMatrix.logAlignPostProb (newSiblingPath);

  const TokSeq newNodeSeq = sibMatrix.sampleParent (newSiblingPath, generator);
  const LogProb logPostNewSeq = sibMatrix.logParentPostProb (newNodeSeq, newSiblingPath);

  const AlignPath oldSiblingPath = Sampler::treePathToSiblingPath (oldAlign.path, history.tree, node);
  const LogProb logPostOldSiblingPath = sibMatrix.logAlignPostProb (oldSiblingPath);

  const TokSeq& oldNodeSeq = validTokenize (oldAlign.ungapped[node].seq, sampler.model.alphabet);
  const LogProb logPostOldSeq = sibMatrix.logParentPostProb (newNodeSeq, oldSiblingPath);

  const AlignPath newNodePath = Sampler::siblingPathToTreePath (newSiblingPath, history.tree, node);

  logForwardProposal = logPostNewSiblingPath + logPostNewSeq;
  logReverseProposal = logPostOldSiblingPath + logPostOldSeq;

  vguard<AlignPath> mergeComponents = { lCladePath, rCladePath, newNodePath };
  AlignPath newPath = alignPathMerge (mergeComponents);

  vguard<FastSeq> newUngapped = oldAlign.ungapped;
  newUngapped[node].seq = detokenize (newNodeSeq, sampler.model.alphabet);

  if (parent >= 0) {  // don't attempt to align to parent if node is root
    const TokSeq pSeq = sampler.removeGapsAndTokenize (history.gapped[parent]);
    const TreeBranchLength pDist = history.tree.branchLength(parent,node);

    const TreeNodeIndex nodeClosestLeaf = history.tree.closestLeaf (node, parent);
    const TreeNodeIndex parentClosestLeaf = history.tree.closestLeaf (parent, node);

    const GuideAlignmentEnvelope branchEnv (sampler.guide.path, parentClosestLeaf, nodeClosestLeaf, sampler.maxDistanceFromGuide);
    const AlignPath& nodeSubtreePath = newPath;
    const vguard<SeqIdx> newNodeEnvPos = guideSeqPos (nodeSubtreePath, node, nodeClosestLeaf);
    const vguard<SeqIdx> oldNodeEnvPos = guideSeqPos (oldAlign.path, node, nodeClosestLeaf);

    const AlignPath pCladePath = Sampler::cladePath (oldAlign.path, history.tree, parent, node);
    const vguard<SeqIdx> parentEnvPos = guideSeqPos (pCladePath, parent, parentClosestLeaf);

    const BranchMatrix newBranchMatrix (sampler.model, pSeq, newNodeSeq, pDist, branchEnv, parentEnvPos, newNodeEnvPos);

    const AlignPath newBranchPath = newBranchMatrix.sample (generator);
    const LogProb logPostNewBranchPath = newBranchMatrix.logPostProb (newBranchPath);

    const BranchMatrix oldBranchMatrix (sampler.model, pSeq, oldNodeSeq, pDist, branchEnv, parentEnvPos, oldNodeEnvPos);
    const AlignPath oldBranchPath = Sampler::treePathToBranchPath (oldAlign.path, history.tree, node);
    const LogProb logPostOldBranchPath = oldBranchMatrix.logPostProb (oldBranchPath);

    logForwardProposal += logPostNewBranchPath;
    logReverseProposal += logPostOldBranchPath;

    mergeComponents.push_back (pCladePath);
    mergeComponents.push_back (newBranchPath);
    newPath = alignPathMerge (mergeComponents);
  }
  
  const Alignment newAlign (newUngapped, newPath);

  newHistory.tree = oldHistory.tree;
  newHistory.gapped = newAlign.gapped();

  oldLogLikelihood = sampler.logLikelihood (oldHistory);
  newLogLikelihood = sampler.logLikelihood (newHistory);
}

Sampler::BranchMatrix::BranchMatrix (const RateModel& model, const TokSeq& xSeq, const TokSeq& ySeq, TreeBranchLength dist, const GuideAlignmentEnvelope& env, const vguard<SeqIdx>& xEnvPos, const vguard<SeqIdx>& yEnvPos)
  : SparseDPMatrix (xSeq, ySeq, env, xEnvPos, yEnvPos),
    model (model),
    probModel (model, dist)
{
  // WRITE ME
}

AlignPath Sampler::BranchMatrix::sample (random_engine& generator) const
{
  // WRITE ME
  return AlignPath();
}

LogProb Sampler::BranchMatrix::logPostProb (const AlignPath& path) const
{
  // WRITE ME
  return -numeric_limits<double>::infinity();
}

Sampler::SiblingMatrix::SiblingMatrix (const RateModel& model, const TokSeq& lSeq, const TokSeq& rSeq, TreeBranchLength plDist, TreeBranchLength prDist, const GuideAlignmentEnvelope& env, const vguard<SeqIdx>& lEnvPos, const vguard<SeqIdx>& rEnvPos)
  : SparseDPMatrix (lSeq, rSeq, env, lEnvPos, rEnvPos),
    model (model),
    lProbModel (model, plDist),
    rProbModel (model, prDist)
{
  // WRITE ME
}

AlignPath Sampler::SiblingMatrix::sampleAlign (random_engine& generator) const
{
  // WRITE ME
  return AlignPath();
}

LogProb Sampler::SiblingMatrix::logAlignPostProb (const AlignPath& plrPath) const
{
  // WRITE ME
  return -numeric_limits<double>::infinity();
}

TokSeq Sampler::SiblingMatrix::sampleParent (const AlignPath& plrPath, random_engine& generator) const
{
  // WRITE ME
  return TokSeq();
}

LogProb Sampler::SiblingMatrix::logParentPostProb (const TokSeq& pSeq, const AlignPath& plrPath) const
{
  // WRITE ME
  return -numeric_limits<double>::infinity();
}

LogProb Sampler::logLikelihood (const History& history) const {
  // WRITE ME
  return -numeric_limits<double>::infinity();
}
