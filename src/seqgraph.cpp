#include <algorithm>
#include "seqgraph.h"
#include "util.h"

SeqGraph::SeqGraph (const Profile& prof, const string& alphabet, const vguard<LogProb>& logInsProb, double minPostProb) {
  const LogProb minLogPostProb = log (minPostProb);
  vguard<vguard<NodeIndex> > stateNodes (prof.size());
  vguard<LogProb> lp (logInsProb.size());
  for (ProfileStateIndex s = 0; s < prof.size(); ++s)
    if (prof.state[s].isNull()) {
      stateNodes[s].push_back (node.size());
      node.push_back (Node());
    } else {
      const LogProb lpNorm = logInnerProduct (logInsProb, prof.state[s].lpAbsorb);
      AlphTok iMax = 0;
      for (AlphTok i = 0; i < logInsProb.size(); ++i) {
	lp[i] = logInsProb[i] + prof.state[s].lpAbsorb[i] - lpNorm;
	if (i == 0 || lp[i] > lp[iMax])
	  iMax = i;
      }
      for (AlphTok i = 0; i < logInsProb.size(); ++i)
	if (i == iMax || lp[i] > minLogPostProb) {
	  stateNodes[s].push_back (node.size());
	  node.push_back (Node());
	  node.back().seq += alphabet[i];
	}
    }
  for (const auto& trans : prof.trans)
    for (auto s : stateNodes[trans.src])
      for (auto d : stateNodes[trans.dest])
	edge.insert (Edge (s, d));
  buildIndices();
  assertToposort();
}

void SeqGraph::buildIndices() {
  for (const auto& e : edge) {
    node[e.src].out.push_back (e);
    node[e.dest].in.push_back (e);
  }
}


void SeqGraph::writeDot (ostream& out) const {
  out << "digraph profile {" << endl;
  for (NodeIndex n = 0; n < nodes(); ++n)
    out << "  n" << n+1 << " [ label = \"" << node[n].seq << "\" ];" << endl;
  for (auto& e : edge)
    out << "  n" << e.src+1 << " -> n" << e.dest+1 << ";" << endl;
  out << "}" << endl;
}

void SeqGraph::assertToposort() const {
  for (auto& e : edge)
    Assert (e.dest > e.src, "SeqGraph is not topologically sorted");
}

vguard<SeqGraph::NodeIndex> SeqGraph::nodeIndices() const {
  vguard<NodeIndex> nl (nodes());
  iota (nl.begin(), nl.end(), 0);
  return nl;
}

vguard<SeqGraph::NodeIndex> SeqGraph::reverseNodeIndices() const {
  auto nl = nodeIndices();
  reverse (nl.begin(), nl.end());
  return nl;
}

SeqGraph SeqGraph::eliminateNull() const {
  map<NodeIndex,set<Edge> > elim;
  set<Edge> keep;
  for (auto src : reverseNodeIndices()) {
    set<Edge> srcOut;
    for (auto& e : node[src].out)
      if (elim.count (e.dest))
	for (auto& e2 : elim[e.dest])
	  srcOut.insert (Edge (src, e2.dest));
      else
	srcOut.insert (e);
    if (node[src].seq.empty())
      elim[src] = srcOut;
    else
      keep.insert (srcOut.begin(), srcOut.end());
  }
  SeqGraph g;
  map<NodeIndex,NodeIndex> old2new;
  for (auto n : nodeIndices())
    if (!node[n].seq.empty()) {
      old2new[n] = g.node.size();
      g.node.push_back (Node());
      g.node.back().seq = node[n].seq;
    }
  for (const auto& e : keep)
    g.edge.insert (Edge (old2new.at(e.src), old2new.at(e.dest)));
  g.buildIndices();
  g.assertToposort();
  return g;
}

SeqGraph SeqGraph::eliminateRedundant() const {
  SeqGraph g;
  return g;
}

SeqGraph SeqGraph::collapseChains() const {
  SeqGraph g;
  return g;
}
