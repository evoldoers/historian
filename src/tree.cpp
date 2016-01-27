#include "tree.h"
#include "knhx.h"
#include "logger.h"

Tree::Tree (const string& nhx) {
  parse (nhx);
}

void Tree::parse (const string& nhx) {
  LogThisAt(5,"Parsing tree: " << nhx << endl);
  ktree_t *tree = kn_parse (nhx.c_str());
  LogThisAt(6,"Tree has " << plural(tree->n,"node") << endl);
  node = vguard<TreeNode> (tree->n);
  for (int n = 0; n < tree->n; ++n) {
    node[n].parent = tree->node[n].parent;
    node[n].child = vector<TreeNodeIndex> (tree->node[n].n);
    for (int c = 0; c < tree->node[n].n; ++c)
      node[n].child[c] = tree->node[n].child[c];
    node[n].name = tree->node[n].name;
    node[n].d = tree->node[n].d;
  }
  kn_free (tree);
}

string Tree::toString (TreeNodeIndex root) const {
  if (isLeaf(root))
    return nodeName(root);
  string s = "(";
  for (int c = 0; c < nChildren(root); ++c) {
    if (c > 0)
      s += ",";
    s += toString(getChild(root,c)) + ":" + to_string(branchLength(getChild(root,c)));
  }
  s += ")";
  return s;
}

string Tree::toString() const {
  return toString (root());
}

string Tree::nodeName (TreeNodeIndex n) const {
  return node[n].name;
}

double Tree::branchLength (TreeNodeIndex n) const {
  return node[n].d;
}

TreeNodeIndex Tree::nodes() const {
  return (TreeNodeIndex) node.size();
}

TreeNodeIndex Tree::root() const {
  return nodes() - 1;
}

TreeNodeIndex Tree::parentNode (TreeNodeIndex n) const {
  return node[n].parent;
}

bool Tree::isLeaf (TreeNodeIndex n) const {
  return node[n].child.empty();
}

size_t Tree::nChildren (TreeNodeIndex n) const {
  return node[n].child.size();
}

TreeNodeIndex Tree::getChild (TreeNodeIndex n, size_t childNum) const {
  return node[n].child[childNum];
}
