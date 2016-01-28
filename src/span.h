#ifndef SPAN_INCLUDED
#define SPAN_INCLUDED

#include <queue>
#include <list>
#include "quickalign.h"
#include "forward.h"

struct AlignGraph {
  struct Edge {
    AlignRowIndex row1, row2;
    AlignRowPath path1, path2;
    LogProb lp;
    Edge() : lp (-numeric_limits<double>::infinity()) { }
    bool operator< (const Edge& e) const { return lp < e.lp; }
    AlignPath path() const;
  };

  struct Partition {
    vguard<size_t> seqSetIdx;
    vguard<set<size_t> > seqSet;
    size_t nSets;
    const set<size_t> empty;
    Partition (size_t n);
    bool inSameSet (const Edge& e) const;
    void merge (const Edge& e);
  };
  
  const vguard<FastSeq>& seqs;
  const RateModel& model;
  const double time;

  vguard<priority_queue<Edge> > edges;
  
  AlignGraph (const vguard<FastSeq>& seqs, const RateModel& model, const double time, ForwardMatrix::random_engine& generator);
  list<AlignPath> minSpanTree();
  AlignPath mstPath();
};

#endif /* SPAN_INCLUDED */

