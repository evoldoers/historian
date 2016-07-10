#ifndef REFINER_INCLUDED
#define REFINER_INCLUDED

#include "sampler.h"

struct Refiner : TreeAlignFuncs {
  typedef DPMatrix::random_engine random_engine;
  
  // Refiner::BranchMatrix
  class BranchMatrix : public BranchMatrixBase {
  public:
    BranchMatrix (const RateModel& model, const PosWeightMatrix& xSeq, const PosWeightMatrix& ySeq, TreeBranchLength dist, const GuideAlignmentEnvelope& env, const vguard<SeqIdx>& xEnvPos, const vguard<SeqIdx>& yEnvPos, AlignRowIndex xRow, AlignRowIndex yRow);

    AlignPath best() const;
  };

  // Refiner member variables
  const RateModel& model;
  int maxDistanceFromGuide;
  
  // Refiner constructor
  Refiner (const RateModel& model);

  // Refiner sampling methods
  inline LogProb logLikelihood (const History& history) const {
    return TreeAlignFuncs::logLikelihood (model, history);
  }
  
  History refine (const History& oldHistory, TreeNodeIndex node) const;
  History refine (const History& oldHistory) const;

  GuideAlignmentEnvelope makeGuide (const Tree& tree, const AlignPath& path, TreeNodeIndex node1, TreeNodeIndex node2) const;
  vguard<SeqIdx> guideSeqPos (const AlignPath& path, AlignRowIndex row) const;
};


#endif /* REFINER_INCLUDED */
