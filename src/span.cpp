#include "span.h"

AlignSpan::AlignSpan (const vguard<FastSeq>& seqs, const RateModel& model, const double time)
  : seqs (seqs),
    model (model),
    time (time)
{
}
