#include "seqgraph.h"

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
      for (auto d : stateNodes[trans.dest]) {
	node[s].out.push_back (edge.size());
	node[d].in.push_back (edge.size());
	edge.push_back (Edge (s, d));
      }
}

void SeqGraph::writeDot (ostream& out) const {
  out << "digraph profile {" << endl;
  for (NodeIndex n = 0; n < nodes(); ++n)
    out << "  n" << n+1 << " [ label = \"" << node[n].seq << "\" ];" << endl;
  for (EdgeIndex e = 0; e < edges(); ++e)
    out << "  n" << edge[e].src+1 << " -> n" << edge[e].dest+1 << ";" << endl;
  out << "}" << endl;
}
