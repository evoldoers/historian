#ifndef SPAN_INCLUDED
#define SPAN_INCLUDED

#include <queue>
#include <list>
#include "quickalign.h"
#include "forward.h"

struct AlignGraph {
  struct TrialEdge {
    AlignRowIndex row1, row2;
    TrialEdge() { }
    TrialEdge (AlignRowIndex src, AlignRowIndex dest)
      : row1(src), row2(dest)
    { }
  };

  struct Edge : TrialEdge {
    LogProb lp;
    Edge() : lp (-numeric_limits<double>::infinity()) { }
    bool operator< (const Edge& e) const { return lp < e.lp; }
  };

  struct Partition {
    vguard<size_t> seqSetIdx;
    vguard<set<size_t> > seqSet;
    size_t nSets;
    Partition (size_t n);
    bool inSameSet (const TrialEdge& e) const;
    void merge (const TrialEdge& e);
  };
  
  const vguard<FastSeq>& seqs;
  const RateModel& model;
  const double time;
  const DiagEnvParams& diagEnvParams;

  vguard<priority_queue<Edge> > edges;
  vguard<map<AlignRowIndex,AlignPath> > edgePath;
  
  AlignGraph (const vguard<FastSeq>& seqs, const RateModel& model, const double time, const DiagEnvParams& diagEnvParams, ForwardMatrix::random_engine& generator);
  AlignGraph (const vguard<FastSeq>& seqs, const RateModel& model, const double time, const DiagEnvParams& diagEnvParams);

  void buildSparseRandomGraph (ForwardMatrix::random_engine& generator);
  void buildDenseGraph();
  void buildGraph (const list<TrialEdge>& trialEdges);

  list<AlignPath> minSpanTree();
  AlignPath mstPath();
  Alignment mstAlign();
  vguard<FastSeq> mstGapped();
};

#endif /* SPAN_INCLUDED */

