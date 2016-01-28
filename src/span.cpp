#include "span.h"

AlignSpan::AlignSpan (const vguard<FastSeq>& seqs, const RateModel& model, const double time, ForwardMatrix::random_engine& generator)
  : seqs (seqs),
    model (model),
    time (time),
    edge (seqs.size())
{
  vguard<size_t> seqSetIdx (seqs.size());
  vguard<set<size_t> > seqSet (seqs.size());
  for (size_t n = 0; n < seqs.size(); ++n) {
    seqSetIdx[n] = n;
    seqSet[n].insert (n);
  }
  size_t nSets = seqs.size();
  
  const size_t nEdges = min ((size_t) (seqs.size() * (seqs.size() - 1) / 2),
			     (size_t) ceil (log(seqs.size()) * (double) seqs.size() / log(2)));
  
  uniform_int_distribution<size_t> dist (0, seqs.size() - 1);
  for (size_t n = 0; n < nEdges || nSets > 1; ++n) {
    size_t src, dest;
    do {
      src = dist (generator);
      dest = dist (generator);
      if (dest < src)
	swap (src, dest);
    } while (src == dest || edge[src].find(dest) != edge[src].end());

    DiagonalEnvelope env (seqs[src], seqs[dest]);
    env.initFull();

    QuickAlignMatrix mx (env, model, time);
    AlignPath path = mx.alignment();

    Edge e;
    e.row1 = path[0];
    e.row2 = path[1];
    e.lp = mx.end;

    edge[src][dest] = e;

    if (seqSetIdx[src] != seqSetIdx[dest]) {
      set<size_t>& destSet = seqSet[seqSetIdx[dest]];
      seqSet[seqSetIdx[src]].insert (destSet.begin(), destSet.end());
      for (auto i : destSet)
	seqSetIdx[i] = seqSetIdx[src];
      destSet.clear();
      --nSets;
    }
  }
}
