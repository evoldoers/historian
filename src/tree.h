#ifndef TREE_INCLUDED
#define TREE_INCLUDED

#include <string>
#include "vguard.h"
#include "fastseq.h"

using namespace std;

#define DefaultNodeNamePrefix "node"

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

  // accessors
  string nodeName (TreeNodeIndex node) const;  // returns empty string if node is unnamed
  TreeBranchLength branchLength (TreeNodeIndex node) const;
  TreeNodeIndex nodes() const;
  TreeNodeIndex root() const;
  TreeNodeIndex parentNode (TreeNodeIndex node) const;
  bool isLeaf (TreeNodeIndex node) const;
  size_t nChildren (TreeNodeIndex node) const;
  TreeNodeIndex getChild (TreeNodeIndex node, size_t childNum) const;
  TreeNodeIndex getSibling (TreeNodeIndex node) const;

  // I/O
  void parse (const string& nhx);
  void validateBranchLengths() const;

  string nodeToString (TreeNodeIndex n) const;  // Newick format, without trailing ";"
  string toString (TreeNodeIndex root) const;  // Newick format, with trailing ";"
  string toString() const;

  void assignInternalNodeNames (const char* prefix = DefaultNodeNamePrefix);
  
  // neighbor-joining
  void buildByNeighborJoining (const vguard<string>& nodeName, const vguard<vguard<TreeBranchLength> >& distanceMatrix);
  void buildByNeighborJoining (const vguard<FastSeq>& seq, const vguard<vguard<TreeBranchLength> >& distanceMatrix);

  // mapping to sequence dataset
  string seqName (TreeNodeIndex node) const;  // guaranteed nonempty
  static string pairParentName (const string& lChildName, double lTime, const string& rChildName, double rTime);

  void reorder (vguard<FastSeq>& seq) const;  // reorders seq so that seq[n].name == seqName(n)
};

#endif /* TREE_INCLUDED */
