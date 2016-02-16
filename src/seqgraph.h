#ifndef SEQGRAPH_INCLUDED
#define SEQGRAPH_INCLUDED

#include "profile.h"

struct SeqGraph {
  typedef size_t NodeIndex;
  typedef size_t EdgeIndex;

  struct Edge {
    NodeIndex src, dest;
    Edge (NodeIndex s, NodeIndex d) : src(s), dest(d) { }
  };

  struct Node {
    vguard<EdgeIndex> in, out;
    string seq;
  };

  vguard<Node> node;
  vguard<Edge> edge;

  SeqGraph (const Profile& prof, const string& alphabet, const vguard<LogProb>& logInsProb, double minPostProb);

  NodeIndex nodes() const { return node.size(); }
  NodeIndex edges() const { return edge.size(); }
  
  void writeDot (ostream& out) const;
};

#endif /* SEQGRAPH_INCLUDED */
