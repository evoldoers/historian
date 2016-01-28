#ifndef SPAN_INCLUDED
#define SPAN_INCLUDED

#include "quickalign.h"
#include "forward.h"

struct AlignSpan {
  struct Edge {
    vguard<bool> row1, row2;
    LogProb lp;
  };

  const vguard<FastSeq>& seqs;
  const RateModel& model;
  const double time;
  vguard<map<AlignRowIndex,Edge> > edge;

  AlignSpan (const vguard<FastSeq>& seqs, const RateModel& model, const double time, ForwardMatrix::random_engine& generator);
  AlignPath minSpanTree() const;
};

#endif /* SPAN_INCLUDED */

