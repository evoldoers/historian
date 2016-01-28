#ifndef SPAN_INCLUDED
#define SPAN_INCLUDED

#include <queue>
#include <list>
#include "quickalign.h"
#include "forward.h"

struct AlignGraph {
  struct Edge {
    AlignRowIndex row1, row2;
    LogProb lp;
    Edge() : lp (-numeric_limits<double>::infinity()) { }
    bool operator< (const Edge& e) const { return lp < e.lp; }
  };

  struct Partition {
    vguard<size_t> seqSetIdx;
    vguard<set<size_t> > seqSet;
    size_t nSets;
    Partition (size_t n);
    bool inSameSet (const Edge& e) const;
    void merge (const Edge& e);
  };
  
  const vguard<FastSeq>& seqs;
  const RateModel& model;
  const double time;

  vguard<priority_queue<Edge> > edges;
  vguard<map<AlignRowIndex,AlignPath> > edgePath;
  
  AlignGraph (const vguard<FastSeq>& seqs, const RateModel& model, const double time, ForwardMatrix::random_engine& generator);
  list<AlignPath> minSpanTree();
  AlignPath mstPath();
};

#endif /* SPAN_INCLUDED */

