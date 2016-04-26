#ifndef SAMPLER_INCLUDED
#define SAMPLER_INCLUDED

#include "model.h"
#include "tree.h"
#include "fastseq.h"

struct SimpleTreePrior {
  double coalescenceRate;
  LogProb treeLogLikelihood (const Tree& tree) const;
};

struct Sampler {
  struct State {
    vguard<FastSeq> gapped;
    Tree tree;
  };

  struct Move {
    enum Type { Branch, TripleBranch, PruneRegraft, NodeHeight };
    Type type;
    TreeNodeIndex node, parent, leftChild, rightChild, oldSibling, newSibling;
    State oldState, newState;
    LogProb logProposal, logInverseProposal, logOldState, logNewState, logHastingsRatio;
  };

  struct AlignmentProposal {
    const RateModel* model;
    FastSeq seq0, seq1;
    TreeBranchLength dist;
    AlignPath path;
    LogProb logPostProb;  // log P(path|model,seq0,seq1,dist)

    AlignmentProposal (const RateModel* model, const FastSeq& seq0, const FastSeq& seq1, TreeBranchLength dist);
    AlignmentProposal (const AlignPath& path, const RateModel* model, const FastSeq& seq0, const FastSeq& seq1, TreeBranchLength dist);
  };

  struct ParentProposal {
    const RateModel* model;
    TreeNodeIndex lNode, rNode, pNode;
    FastSeq lSeq, rSeq, pSeq;
    TreeBranchLength plDist, prDist;
    AlignPath plrPath;
    LogProb logPostProb;  // log P(pSeq,plrPath|model,lSeq,rSeq,plDist,prDist,lrPath)

    ParentProposal (const RateModel* model, const FastSeq& lSeq, const FastSeq& rSeq, TreeBranchLength plDist, TreeBranchLength prDist, const AlignPath& lrPath);
    ParentProposal (const FastSeq& pSeq, const AlignPath& plrPath, const RateModel* model, const FastSeq& lSeq, const FastSeq& rSeq, TreeBranchLength plDist, TreeBranchLength prDist, const AlignPath& lrPath);
};
  
  RateModel model;
  SimpleTreePrior treePrior;
};

#endif /* SAMPLER_INCLUDED */
