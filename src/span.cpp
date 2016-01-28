#include "span.h"

AlignGraph::Partition::Partition (size_t n)
  : seqSetIdx (n),
    seqSet (n),
    nSets (n)
{
  for (size_t i = 0; i < n; ++i) {
    seqSetIdx[i] = i;
    seqSet[i].insert (i);
  }
}

bool AlignGraph::Partition::inSameSet (const AlignGraph::Edge& e) const {
  return seqSetIdx[e.row1] == seqSetIdx[e.row2];
}

void AlignGraph::Partition::merge (const AlignGraph::Edge& e) {
  if (!inSameSet(e)) {
    const size_t idx1 = seqSetIdx[e.row1];
    const size_t idx2 = seqSetIdx[e.row2];
    set<size_t>& set1 = seqSet[idx1];
    set<size_t>& set2 = seqSet[idx2];
    for (auto n2 : set2)
      seqSetIdx[n2] = idx1;
    set1.insert (set2.begin(), set2.end());
    set2.clear();
    --nSets;
  }
}

AlignPath AlignGraph::Edge::path() const {
  AlignPath p;
  p[row1] = path1;
  p[row2] = path2;
  return p;
}

AlignGraph::AlignGraph (const vguard<FastSeq>& seqs, const RateModel& model, const double time, ForwardMatrix::random_engine& generator)
  : seqs (seqs),
    model (model),
    time (time),
    edges (seqs.size())
{
  Partition part (seqs.size());
  vguard<set<AlignRowIndex> > edgeDests (seqs.size());

  const size_t nEdges = min ((size_t) (seqs.size() * (seqs.size() - 1) / 2),
			     (size_t) ceil (log(seqs.size()) * (double) seqs.size() / log(2)));
  
  uniform_int_distribution<size_t> dist (0, seqs.size() - 1);
  for (size_t n = 0; n < nEdges || part.nSets > 1; ++n) {
    size_t src, dest;
    do {
      src = dist (generator);
      dest = dist (generator);
      if (dest < src)
	swap (src, dest);
    } while (src == dest || edgeDests[src].count(dest));

    DiagonalEnvelope env (seqs[src], seqs[dest]);
    env.initFull();

    QuickAlignMatrix mx (env, model, time);
    AlignPath path = mx.alignment();

    Edge e;
    e.row1 = src;
    e.row2 = dest;
    e.path1 = path[0];
    e.path2 = path[1];
    e.lp = mx.end;

    edgeDests[src].insert (dest);
    edges[src].push (e);

    part.merge (e);
  }
}

list<AlignPath> AlignGraph::minSpanTree() {
  list<AlignPath> paths;
  Partition part (seqs.size());
  while (part.nSets > 1) {
    Edge best;
    for (auto src : part.seqSet.front()) {
      while (!edges[src].empty() && part.inSameSet (edges[src].top()))
	edges[src].pop();
      if (!edges[src].empty() && edges[src].top().lp > best.lp)
	best = edges[src].top();
    }
    paths.push_back (best.path());
    part.merge (best);
  }
  edges.clear();  // this method is unambiguously destructive...
  return paths;
}

AlignPath AlignGraph::mstPath() {
  const list<AlignPath> pathList = minSpanTree();
  const vguard<AlignPath> pathVec (pathList.begin(), pathList.end());
  return alignPathMerge (pathVec);
}
