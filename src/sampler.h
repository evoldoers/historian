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
    State oldState, newState;
    LogProb logHastingsRatio;
  };
  
  RateModel model;
  SimpleTreePrior treePrior;

  Move realignBranch (const State& oldState, TreeNodeIndex childNode) const;
  Move realignNodeNeighborhood (const State& oldState, TreeNodeIndex node) const;
  Move pruneAndRegraft (const State& oldState, TreeNodeIndex childNode, TreeNodeIndex newSibling) const;
  Move resampleNodeHeight (const State& oldState, TreeNodeIndex node) const;
};

#endif /* SAMPLER_INCLUDED */
