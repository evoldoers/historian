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

  typedef DPMatrix::random_engine random_engine;

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

  struct AlignmentMatrix {
    const RateModel* model;
    FastSeq parentSeq, childSeq;
    TreeBranchLength dist;

    AlignmentMatrix (const RateModel* model, const FastSeq& seq0, const FastSeq& seq1, TreeBranchLength dist);
    void sample (AlignPath& path, random_engine& generator) const;
    LogProb logPostProb (const AlignPath& path) const;

    static inline AlignRowIndex parent() const { return 0; }
    static inline AlignRowIndex child() const { return 1; }
  };

  struct ParentMatrix {
    const RateModel* model;
    FastSeq lSeq, rSeq;
    TreeBranchLength plDist, prDist;
    AlignPath lrPath;

    ParentMatrix (const RateModel* model, const FastSeq& lSeq, const FastSeq& rSeq, TreeBranchLength plDist, TreeBranchLength prDist, const AlignPath& lrPath);
    void sample (FastSeq& pSeq, AlignPath& plrPath, random_engine& generator) const;
    LogProb logPostProb (const FastSeq& pSeq, const AlignPath& plrPath) const;

    static inline AlignRowIndex lChild() const { return 0; }
    static inline AlignRowIndex rChild() const { return 1; }
    static inline AlignRowIndex parent() const { return 2; }
  };
  
  RateModel model;
  SimpleTreePrior treePrior;
};

#endif /* SAMPLER_INCLUDED */
