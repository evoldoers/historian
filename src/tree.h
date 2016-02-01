#ifndef TREE_INCLUDED
#define TREE_INCLUDED

#include <string>
#include "vguard.h"
#include "fastseq.h"

using namespace std;

typedef int TreeNodeIndex;
typedef double TreeBranchLength;

struct TreeNode {
  TreeNodeIndex parent;
  vguard<TreeNodeIndex> child;
  string name;
  TreeBranchLength d;
};

struct Tree {
  vector<TreeNode> node;

  Tree() { }
  Tree (const string& nhx);

  string nodeName (TreeNodeIndex node) const;
  TreeBranchLength branchLength (TreeNodeIndex node) const;
  TreeNodeIndex nodes() const;
  TreeNodeIndex root() const;
  TreeNodeIndex parentNode (TreeNodeIndex node) const;
  bool isLeaf (TreeNodeIndex node) const;
  size_t nChildren (TreeNodeIndex node) const;
  TreeNodeIndex getChild (TreeNodeIndex node, size_t childNum) const;
  TreeNodeIndex getSibling (TreeNodeIndex node) const;

  void parse (const string& nhx);
  
  string nodeToString (TreeNodeIndex n) const;  // omits trailing ";"
  string toString (TreeNodeIndex root) const;
  string toString() const;

  void buildByNeighborJoining (const vguard<string>& nodeName, const vguard<vguard<TreeBranchLength> >& distanceMatrix);
  void buildByNeighborJoining (const vguard<FastSeq>& seq, const vguard<vguard<TreeBranchLength> >& distanceMatrix);
};

#endif /* TREE_INCLUDED */
