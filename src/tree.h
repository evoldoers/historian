#ifndef TREE_INCLUDED
#define TREE_INCLUDED

#include <string>
#include <set>
#include "vguard.h"
#include "fastseq.h"

using namespace std;

#define DefaultNodeNamePrefix "node"
#define DefaultNewRootName "root"

typedef int TreeNodeIndex;
typedef double TreeBranchLength;

struct TreeNode {
  TreeNodeIndex parent;
  vguard<TreeNodeIndex> child;
  string name;
  TreeBranchLength d;
};

#define TREE_MIN_BRANCH_LEN 1e-6
struct Tree {
  vector<TreeNode> node;

  static double minBranchLength;
  
  Tree() { }
  Tree (const string& nhx);

  // accessors
  string nodeName (TreeNodeIndex node) const;  // returns empty string if node is unnamed
  TreeBranchLength branchLength (TreeNodeIndex node) const;
  TreeBranchLength branchLength (TreeNodeIndex node1, TreeNodeIndex node2) const;
  TreeNodeIndex nodes() const;
  TreeNodeIndex root() const;
  TreeNodeIndex parentNode (TreeNodeIndex node) const;
  bool isLeaf (TreeNodeIndex node) const;
  size_t nChildren (TreeNodeIndex node) const;
  TreeNodeIndex getChild (TreeNodeIndex node, size_t childNum) const;
  TreeNodeIndex getSibling (TreeNodeIndex node) const;
  vguard<TreeNodeIndex> getSiblings (TreeNodeIndex node) const;

  set<TreeNodeIndex> nodeAndAncestors (TreeNodeIndex node) const;
  set<TreeNodeIndex> nodeAndDescendants (TreeNodeIndex node) const;
  TreeNodeIndex mostRecentCommonAncestor (TreeNodeIndex node1, TreeNodeIndex node2) const;
  
  TreeNodeIndex findNode (const string& name) const;

  bool isBinary() const;
  void assertBinary() const;

  bool isUltrametric (double epsilon = 1e-4) const;
  void assertUltrametric (double epsilon = 1e-4) const;

  bool isPostorderSorted() const;
  void assertPostorderSorted() const;

  vguard<TreeNodeIndex> rerootedChildren (TreeNodeIndex node, TreeNodeIndex parent) const;
  vguard<TreeNodeIndex> rerootedPreorderSort (TreeNodeIndex newRoot, TreeNodeIndex parentOfRoot = -1) const;
  vguard<TreeNodeIndex> rerootedParent (TreeNodeIndex newRoot) const;

  vguard<TreeNodeIndex> preorderSort() const;
  vguard<TreeNodeIndex> postorderSort() const;

  Tree reorderNodes (const vguard<TreeNodeIndex>& newOrder) const;
  void detach (TreeNodeIndex node);
  void setParent (TreeNodeIndex node, TreeNodeIndex parent, TreeBranchLength branchLength);  // WARNING! does not check for cycles, may leave tree in a non-preorder-sorted state

  vguard<TreeBranchLength> distanceFromRoot() const;
  
  // I/O
  void parse (const string& nhx);
  void validateBranchLengths() const;

  static string branchLengthString (TreeBranchLength d);

  pair<string,TreeBranchLength> nodeDescriptor (TreeNodeIndex n, TreeNodeIndex parent) const;
  string nodeToString (TreeNodeIndex n, TreeNodeIndex parent) const;  // Newick format, without trailing ";"
  string nodeToString (TreeNodeIndex n) const;  // does not re-root
  string toString (TreeNodeIndex root, TreeNodeIndex parent) const;  // Newick format, with trailing ";"
  string toString (TreeNodeIndex root) const;  // does not re-root
  string toString() const;

  string toStringRerootedAbove (TreeNodeIndex node, const char* newRootName = DefaultNewRootName) const;
  Tree rerootAbove (TreeNodeIndex node, const char* newRootName = DefaultNewRootName) const;
  Tree rerootAbove (const string& name, const char* newRootName = DefaultNewRootName) const;

  // neighbor-joining & UPGMA algorithms
  void buildByNeighborJoining (const vguard<string>& nodeName, const vguard<vguard<TreeBranchLength> >& distanceMatrix);
  void buildByNeighborJoining (const vguard<FastSeq>& seq, const vguard<vguard<TreeBranchLength> >& distanceMatrix);
  void buildByUPGMA (const vguard<string>& nodeName, const vguard<vguard<TreeBranchLength> >& distanceMatrix);
  void buildByUPGMA (const vguard<FastSeq>& seq, const vguard<vguard<TreeBranchLength> >& distanceMatrix);

  // mapping to sequence dataset
  string seqName (TreeNodeIndex node) const;  // guaranteed nonempty
  static string pairParentName (const string& lChildName, double lTime, const string& rChildName, double rTime);

  bool allNodesNamed() const;  // true if all nodes have nonempty names
  void assertAllNodesNamed() const;
  
  void reorderSeqs (vguard<FastSeq>& seq) const;  // reorders seq so that seq[n].name == seqName(n)

  bool nodesMatchSeqs (const vguard<FastSeq>& seq) const;  // true if seq[n].name == nodeName(n)
  void assertNodesMatchSeqs (const vguard<FastSeq>& seq) const;

  void assignInternalNodeNames (const char* prefix = DefaultNodeNamePrefix);
  void assignInternalNodeNames (vguard<FastSeq>& seq, const char* prefix = DefaultNodeNamePrefix);

  // general helpers
  
  TreeNodeIndex closestLeaf (TreeNodeIndex node, TreeNodeIndex parent = -1) const;
};

#endif /* TREE_INCLUDED */
