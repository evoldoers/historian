#ifndef TREE_INCLUDED
#define TREE_INCLUDED

#include <string>
#include "vguard.h"

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
  double branchLength (TreeNodeIndex node) const;
  TreeNodeIndex nodes() const;
  TreeNodeIndex root() const;
  TreeNodeIndex parentNode (TreeNodeIndex node) const;
  bool isLeaf (TreeNodeIndex node) const;
  size_t nChildren (TreeNodeIndex node) const;
  TreeNodeIndex getChild (TreeNodeIndex node, size_t childNum) const;

  void parse (const string& nhx);
  
  string toString (TreeNodeIndex root) const;
  string toString() const;
};

#endif /* TREE_INCLUDED */
