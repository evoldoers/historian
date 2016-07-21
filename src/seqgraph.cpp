#include <algorithm>
#include "seqgraph.h"
#include "util.h"
#include "logger.h"

SeqGraph::SeqGraph (const Profile& prof, const string& alphabet, const vguard<LogProb>& logCptWeight, const vguard<vguard<LogProb> >& logInsProb, double minPostProb) {
  const LogProb minLogPostProb = log (minPostProb);
  vguard<vguard<NodeIndex> > stateNodes (prof.size());
  for (ProfileStateIndex s = 0; s < prof.size(); ++s)
    if (prof.state[s].isNull()) {
      stateNodes[s].push_back (node.size());
      node.push_back (Node());
    } else {
      LogProb lpNorm = -numeric_limits<double>::infinity();
      for (int cpt = 0; cpt < (int) logInsProb.size(); ++cpt)
	log_accum_exp (lpNorm, logCptWeight[cpt] + logInnerProduct (logInsProb[cpt], prof.state[s].lpAbsorb[cpt]));
      AlphTok iMax = 0;
      vguard<LogProb> lp (alphabet.size(), -numeric_limits<double>::infinity());
      for (AlphTok i = 0; i < (int) lp.size(); ++i) {
	for (int cpt = 0; cpt < (int) logInsProb.size(); ++cpt)
	  log_accum_exp (lp[i], logCptWeight[cpt] + logInsProb[cpt][i] + prof.state[s].lpAbsorb[cpt][i] - lpNorm);
	if (i == 0 || lp[i] > lp[iMax])
	  iMax = i;
      }
      for (AlphTok i = 0; i < alphabet.size(); ++i)
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
}

void SeqGraph::buildIndices() {
  for (const auto& e : edge) {
    node[e.src].out.push_back (e);
    node[e.dest].in.push_back (e);
  }
  LogThisAt(3,"Sequence graph has " << plural(nodes(),"node") << " and " << plural(edges(),"edge") << endl);
  assertToposort();
}

void SeqGraph::writeDot (ostream& out) const {
  out << "digraph profile {" << endl;
  for (NodeIndex n = 0; n < nodes(); ++n)
    out << "  n" << n+1 << " [ shape = rect, label = \"" << node[n].seq << "\" ];" << endl;
  for (auto& e : edge)
    out << "  n" << e.src+1 << " -> n" << e.dest+1 << ";" << endl;
  out << "}" << endl;
}

void SeqGraph::assertToposort() const {
  for (auto& e : edge) {
    Assert (e.dest != e.src, "SeqGraph has a self-looping state");
    Assert (e.dest > e.src, "SeqGraph is not topologically sorted");
  }
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
    if (node[src].isNull())
      elim[src] = srcOut;
    else
      keep.insert (srcOut.begin(), srcOut.end());
  }
  SeqGraph g;
  if (elim.empty())
    g = *this;
  else {
    map<NodeIndex,NodeIndex> old2new;
    for (auto n : nodeIndices())
      if (!node[n].isNull()) {
	old2new[n] = g.node.size();
	g.node.push_back (Node());
	g.node.back().seq = node[n].seq;
      }
    for (const auto& e : keep)
      g.edge.insert (Edge (old2new.at(e.src), old2new.at(e.dest)));
    LogThisAt(3,"Eliminated " << plural(nodes() - g.nodes(),"null node") << endl);
    g.buildIndices();
  }
  return g;
}

SeqGraph SeqGraph::eliminateDuplicates() const {

  struct SeqGraphNodeSummary {
    string seq;
    set<NodeIndex> dest;
    SeqGraphNodeSummary (const Node& node, const map<NodeIndex,NodeIndex>& equiv)
      : seq (node.seq)
    {
      for (const auto& e : node.out)
	dest.insert (equiv.count(e.dest) ? equiv.at(e.dest) : e.dest);
    }
    bool operator== (const SeqGraphNodeSummary& summ) const {
      return seq == summ.seq && dest == summ.dest;
    }
    bool operator< (const SeqGraphNodeSummary& summ) const {
      return seq == summ.seq ? (dest < summ.dest) : (seq < summ.seq);
    }
    string str() const {
      return seq + " " + to_string_join (dest);
    }
  };

  map<NodeIndex,NodeIndex> equiv, old2new;
  map<SeqGraphNodeSummary,NodeIndex> unique;
  for (auto n : reverseNodeIndices()) {
    const SeqGraphNodeSummary summ (node[n], equiv);
    if (unique.count(summ))
      equiv[n] = unique[summ];
    else
      unique[summ] = n;
  }
  SeqGraph g;
  if (equiv.empty())
    g = *this;
  else {
    for (auto n : nodeIndices())
      if (!equiv.count(n)) {
	old2new[n] = g.node.size();
	g.node.push_back (Node());
	g.node.back().seq = node[n].seq;
      }
    for (auto& e : edge)
      if (old2new.count(e.src))
	g.edge.insert (Edge (old2new.at(e.src),
			     old2new.at (equiv.count(e.dest) ? equiv.at(e.dest) : e.dest)));
    LogThisAt(3,"Eliminated " << plural(nodes() - g.nodes(),"duplicate node") << endl);
    g.buildIndices();
  }
  return g;
}

SeqGraph SeqGraph::collapseChains() const {
  map<NodeIndex,NodeIndex> chainEnd, old2new;
  map<NodeIndex,string> chainSeq;
  set<NodeIndex> elim;
  NodeIndex dest;
  for (auto n : reverseNodeIndices())
    if (node[n].out.size() == 1
	&& chainEnd.count(dest = node[n].out.front().dest)
	&& node[dest].in.size() == 1) {
      chainEnd[n] = chainEnd.at(dest);
      chainSeq[chainEnd[n]] = node[n].seq + chainSeq[chainEnd[n]];
      elim.insert (n);
    } else if (node[n].in.size() == 1) {
      chainEnd[n] = n;
      chainSeq[n] = node[n].seq;
    }
  SeqGraph g;
  if (elim.empty())
    g = *this;
  else {
    for (auto n : nodeIndices())
      if (!elim.count(n)) {
	old2new[n] = g.node.size();
	g.node.push_back (Node());
	g.node.back().seq = chainSeq.count(n) ? chainSeq.at(n) : node[n].seq;
      }
    for (auto& e : edge)
      if (old2new.count (e.src))
	g.edge.insert (Edge (old2new.at(e.src),
			     old2new.at (chainEnd.count(e.dest) ? chainEnd.at(e.dest) : e.dest)));
    LogThisAt(3,"Eliminated " << plural(nodes() - g.nodes(),"chained node") << endl);
    g.buildIndices();
  }
  return g;
}

SeqGraph SeqGraph::mergeCharClasses() const {

  struct SeqGraphNodeSummary {
    set<NodeIndex> src, dest;
    SeqGraphNodeSummary (const Node& node, const map<NodeIndex,NodeIndex>& equiv)
    {
      for (const auto& e : node.in)
	src.insert (equiv.count(e.src) ? equiv.at(e.src) : e.src);
      for (const auto& e : node.out)
	dest.insert (equiv.count(e.dest) ? equiv.at(e.dest) : e.dest);
    }
    bool operator== (const SeqGraphNodeSummary& summ) const {
      return src == summ.src && dest == summ.dest;
    }
    bool operator< (const SeqGraphNodeSummary& summ) const {
      return src == summ.src ? (dest < summ.dest) : (src < summ.src);
    }
    string str() const {
      return to_string_join (src) + " " + to_string_join (dest);
    }
  };

  map<NodeIndex,NodeIndex> equiv, old2new;
  map<SeqGraphNodeSummary,NodeIndex> classRep;
  map<NodeIndex,string> classChars;
  for (auto n : reverseNodeIndices())
    if (node[n].seq.size() == 1) {
      const SeqGraphNodeSummary summ (node[n], equiv);
      if (classRep.count(summ)) {
	equiv[n] = classRep[summ];
	classChars[classRep[summ]] = node[n].seq + classChars[classRep[summ]];
      } else {
	classRep[summ] = n;
	classChars[n] = node[n].seq;
      }
    }

  SeqGraph g;
  if (equiv.empty())
    g = *this;
  else {
    for (auto n : nodeIndices())
      if (!equiv.count(n)) {
	old2new[n] = g.node.size();
	g.node.push_back (Node());
	if (classChars.count(n) && classChars[n].size() > 1)
	  (g.node.back().seq += "[") += classChars[n] + "]";
	else
	  g.node.back().seq = node[n].seq;
      }
    for (auto& e : edge)
      if (old2new.count(e.src) && old2new.count(e.dest))
	g.edge.insert (Edge (old2new[e.src], old2new[e.dest]));
    LogThisAt(3,"Eliminated " << plural(nodes() - g.nodes(),"class node") << endl);
    g.buildIndices();
  }

  return g;
}

SeqGraph SeqGraph::simplify() const {
  return eliminateNull().eliminateDuplicates().mergeCharClasses().collapseChains();
}
