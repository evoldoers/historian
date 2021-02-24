#include <cmath>
#include <gsl/gsl_math.h>
#include <ios>

#include "tree.h"
#include "knhx.h"
#include "logger.h"

double Tree::minBranchLength = TREE_MIN_BRANCH_LEN;

Tree::Tree (const string& nhx) {
  parse (nhx);
}

void Tree::parse (const string& nhx) {
  LogThisAt(8,"Parsing tree: " << nhx << endl);
  ktree_t *tree = kn_parse (nhx.c_str());
  LogThisAt(8,"Tree has " << plural(tree->n,"node") << endl);
  node = vguard<TreeNode> (tree->n);
  set<string> names;
  for (int n = 0; n < tree->n; ++n) {
    node[n].parent = tree->node[n].parent;
    node[n].child = vector<TreeNodeIndex> (tree->node[n].n);
    for (int c = 0; c < tree->node[n].n; ++c)
      node[n].child[c] = tree->node[n].child[c];
    node[n].name = tree->node[n].name;
    if (tree->node[n].d >= 0)
      node[n].d = max (tree->node[n].d, minBranchLength);
    else
      node[n].d = tree->node[n].d;
    if (node[n].name.size()) {
      Require (names.count (node[n].name) == 0, "Duplicate node name '%s' in tree: %s", node[n].name.c_str(), nhx.c_str());
      names.insert (node[n].name);
    }
  }
  kn_free (tree);
}

void Tree::validateBranchLengths() const {
  for (size_t n = 0; n + 1 < node.size(); ++n) {
    Require (branchLength(n) >= 0, "Node in tree is missing branch length: %s", seqName(n).c_str());
    Require (branchLength(n) >= minBranchLength, "Node in tree has a lower-than-minimal branch length: %s", seqName(n).c_str());
  }
}

vguard<TreeNodeIndex> Tree::rerootedChildren (TreeNodeIndex node, TreeNodeIndex parent) const {
  vguard<TreeNodeIndex> children;
  for (int c = 0; c < nChildren(node); ++c)
    if (getChild(node,c) != parent)
      children.push_back (getChild(node,c));
  if (parentNode(node) >= 0 && parentNode(node) != parent)
    children.push_back (parentNode(node));
  return children;
}

string Tree::branchLengthString (TreeBranchLength d) {
  ostringstream ds;
  if (d >= 0) {
    ds.unsetf(std::ios_base::floatfield);
    ds << ':' << d;
  }
  return ds.str();
}

pair<string,TreeBranchLength> Tree::nodeDescriptor (TreeNodeIndex n, TreeNodeIndex parent) const {
  const vguard<TreeNodeIndex> children = rerootedChildren (n, parent);
  if (children.empty())
    return pair<string,TreeBranchLength> (nodeName(n), branchLength(parent,n));
  if (children.size() == 1) {
    pair<string,TreeBranchLength> cd = nodeDescriptor (children.front(), n);
    return pair<string,TreeBranchLength> (cd.first, cd.second + branchLength(parent,n));
  }
  string s = "(";
  for (size_t c = 0; c < children.size(); ++c) {
    pair<string,TreeBranchLength> cd = nodeDescriptor (children[c], n);
    if (c > 0)
      s += ",";
    s += cd.first + branchLengthString (cd.second);
  }
  s += ")" + nodeName(n);
  return pair<string,TreeBranchLength> (s, branchLength(parent,n));
}

string Tree::nodeToString (TreeNodeIndex root) const {
  return nodeToString (root, parentNode(root));
}

string Tree::nodeToString (TreeNodeIndex root, TreeNodeIndex parent) const {
  return nodeDescriptor(root,parent).first;
}

string Tree::toString() const {
  return toString(root());
}

string Tree::toString (TreeNodeIndex n, TreeNodeIndex p) const {
  return nodeToString(n,p) + ";";
}

string Tree::toString (TreeNodeIndex n) const {
  return toString (n, parentNode(n));
}

string Tree::toStringRerootedAbove (TreeNodeIndex node, const char* newRootName) const {
  if (node == root() || parentNode(node) == root())
    return toString();
  const TreeNodeIndex parent = parentNode(node);
  auto nd = nodeDescriptor(node,parent), pd = nodeDescriptor(parent,node);
  return string("(") + nd.first + branchLengthString(nd.second/2) + "," + pd.first + branchLengthString(pd.second/2) + ")" + newRootName + ";";
}

Tree Tree::rerootAbove (TreeNodeIndex node, const char* newRootName) const {
  return Tree (toStringRerootedAbove (node, newRootName));
}

Tree Tree::rerootAbove (const string& name, const char* newRootName) const {
  return rerootAbove (findNode (name), newRootName);
}

string Tree::nodeName (TreeNodeIndex n) const {
  return node[n].name;
}

TreeBranchLength Tree::branchLength (TreeNodeIndex n) const {
  return node[n].d;
}

TreeBranchLength Tree::branchLength (TreeNodeIndex node1, TreeNodeIndex node2) const {
  if (node1 == parentNode(node2))
    return branchLength(node2);
  else if (node2 == parentNode(node1))
    return branchLength(node1);
  Abort ("Nodes %d and %d are not connected by a branch", node1, node2);
  return -1;
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

TreeNodeIndex Tree::getSibling (TreeNodeIndex node) const {
  const TreeNodeIndex parent = parentNode(node);
  Assert (parent >= 0, "Attempt to find sibling of root node");
  Assert (nChildren(parent) == 2, "Attempt to find sibling in non-binary tree");
  return getChild(parent,0) == node ? getChild(parent,1) : getChild(parent,0);
}

vguard<TreeNodeIndex> Tree::getSiblings (TreeNodeIndex n) const {
  vguard<TreeNodeIndex> sibs;
  const TreeNodeIndex parent = parentNode(n);
  if (parent >= 0)
    for (auto s: node[parent].child)
      if (s != n)
	sibs.push_back (s);
  return sibs;
}

bool Tree::isBinary() const {
  for (TreeNodeIndex node = 0; node < nodes(); ++node)
    if (!isLeaf(node) && nChildren(node) != 2)
      return false;
  return true;
}

void Tree::assertBinary() const {
  for (TreeNodeIndex node = 0; node < nodes(); ++node)
    if (!isLeaf(node))
      Assert (nChildren(node) == 2, "Tree is not binary: node %d has %s\nSubtree rooted at %d: %s\n", node, plural(nChildren(node),"child","children").c_str(), node, toString(node).c_str());
}

bool Tree::isUltrametric (double epsilon) const {
  const auto dist = distanceFromRoot();
  TreeBranchLength minDist = numeric_limits<double>::infinity();
  for (TreeNodeIndex node = 0; node < nodes(); ++node)
    if (isLeaf(node))
      minDist = min (minDist, dist[node]);
  for (TreeNodeIndex node = 0; node < nodes(); ++node)
    if (isLeaf(node))
      if (gsl_fcmp (dist[node], minDist, epsilon) != 0)
	return false;
  return true;
}

void Tree::assertUltrametric (double epsilon) const {
  if (!isUltrametric(epsilon)) {
    const auto dist = distanceFromRoot();
    for (TreeNodeIndex node = 0; node < nodes(); ++node)
      if (isLeaf(node))
	Warn ("Node #%d (%s) distance from root is %g", node, nodeName(node).c_str(), dist[node]);
    Abort ("Tree is not ultrametric");
  }
}

bool Tree::isPostorderSorted() const {
  for (TreeNodeIndex n = 0; n < root(); ++n)
    if (parentNode(n) <= n)
      return false;
  return parentNode(root()) < 0;
}

void Tree::assertPostorderSorted() const {
  Assert (isPostorderSorted(), "Tree nodes are not sorted in postorder");
}

TreeNodeIndex Tree::findNode (const string& name) const {
  for (TreeNodeIndex n = 0; n < nodes(); ++n)
    if (nodeName(n) == name)
      return n;
  Abort ("Couldn't find tree node %s", name.c_str());
  return 0;
}

bool Tree::hasNode (const string& name) const {
  for (TreeNodeIndex n = 0; n < nodes(); ++n)
    if (nodeName(n) == name)
      return true;
  return false;
}

void Tree::buildByNeighborJoining (const vguard<string>& nodeName, const vguard<vguard<TreeBranchLength> >& distanceMatrix) {
  // check that there are more than 2 nodes
  Assert (nodeName.size() >= 2, "Fewer than 2 nodes; can't make a binary tree");
  // clear the existing tree
  node.clear();
  // copy distance matrix
  vguard<vguard<TreeBranchLength> > dist = distanceMatrix;
  // estimate tree by neighbor-joining
  // algorithm follows description in Durbin et al, pp170-171
  // first, initialise the list of active nodes
  set<TreeNodeIndex> activeNodes;
  for (TreeNodeIndex n = 0; n < (int) nodeName.size(); ++n) {
    activeNodes.insert (n);
    node.push_back (TreeNode());
    node.back().name = nodeName[n];
    node.back().parent = -1;
  }
  // main loop
  vguard<TreeBranchLength> avgDist;
  while (true)
    {
      // get number of active nodes
      const int nActiveNodes = activeNodes.size();
      // loop exit test
      if (nActiveNodes == 2) break;
      Assert (nActiveNodes > 2, "Fewer than 2 nodes left -- should never get here");
      // calculate average distances from each node
      avgDist = vguard<TreeBranchLength> (nodes(), (double) 0);
      for (auto ni : activeNodes)
	{
	  double a_i = 0;
	  for (auto nj : activeNodes)
	    if (nj != ni)
	      a_i += dist[ni][nj];
	  avgDist[ni] = a_i / (double) (nActiveNodes - 2);
	  LogThisAt(7,"Distance correction for node " << ni << " is " << avgDist[ni] << endl);
	}
      // find minimal compensated distance (with avg distances subtracted off)
      bool isFirstPair = true;
      TreeBranchLength minDist = 0;
      TreeNodeIndex min_i = -1, min_j = -1;
      set<TreeNodeIndex>::const_iterator node_i_p, node_j_p, activeNodes_end_p;
      activeNodes_end_p = activeNodes.end();
      for (node_i_p = activeNodes.begin(); node_i_p != activeNodes_end_p; ++node_i_p)
	for (node_j_p = node_i_p, ++node_j_p; node_j_p != activeNodes_end_p; ++node_j_p)
	  {
	    const double compensatedDist = dist [*node_i_p] [*node_j_p] - avgDist [*node_i_p] - avgDist [*node_j_p];
	    LogThisAt(7,"Compensated distance from node " << *node_i_p << " to node " << *node_j_p << " is " << compensatedDist << endl);
	    if (isFirstPair || compensatedDist < minDist)
	      {
		min_i = *node_i_p;
		min_j = *node_j_p;
		minDist = compensatedDist;
		isFirstPair = false;
	      }
	  }
      // nodes min_i and min_j are neighbors -- join them with new index k
      // first, calculate new distances as per NJ algorithm
      const TreeNodeIndex k = nodes();
      dist.push_back (vguard<TreeBranchLength> (k + 1));
      dist[k][k] = 0;
      const TreeBranchLength d_ij = dist[min_i][min_j];
      for (TreeNodeIndex m = 0; m < k; ++m)
	dist[m].push_back (dist[k][m] = 0.5 * (dist[min_i][m] + dist[min_j][m] - d_ij));
      TreeBranchLength d_ik = 0.5 * (d_ij + avgDist[min_i] - avgDist[min_j]);
      TreeBranchLength d_jk = d_ij - d_ik;
      LogThisAt(8,"Before Kuhner-Felsenstein:\ni=" << min_i << ", j=" << min_j << ", k=" << k << ", d_ij=" << d_ij << ", d_ik=" << d_ik << ", d_jk=" << d_jk << "\nDistances from k to other nodes: " << to_string_join(dist[k]) << endl);
      // apply Kuhner-Felsenstein correction to prevent negative branch lengths
      // also enforce minimum branch lengths here
      if (d_ik < minBranchLength)
	{
	  d_jk -= d_ik - minBranchLength;
	  d_ik = minBranchLength;
	}
      if (d_jk < 0)
	{
	  d_ik -= d_jk - minBranchLength;
	  d_jk = minBranchLength;
	}
      dist[min_i][k] = dist[k][min_i] = d_ik;
      dist[min_j][k] = dist[k][min_j] = d_jk;
      // now update the Tree
      node.push_back (TreeNode());
      node[k].child.push_back (min_i);
      node[k].child.push_back (min_j);
      node[min_i].parent = k;
      node[min_i].d = max (0., d_ik);
      node[min_j].parent = k;
      node[min_j].d = max (0., d_jk);
      LogThisAt(7,"Joining nodes " << min_i << " and " << min_j << " to common ancestor " << k << " (branch lengths: " << k << "->" << min_i << " = " << d_ik << ", " << k << "->" << min_j << " = " << d_jk << ")" << endl);
      activeNodes.erase (min_i);
      activeNodes.erase (min_j);
      activeNodes.insert (k);
    }
  // make the root node
  set<TreeNodeIndex>::iterator iter = activeNodes.begin();
  const TreeNodeIndex i = *iter;
  const TreeNodeIndex j = *++iter;
  const double d = max (dist[i][j], 0.);  // don't correct the last node, just keep branch length non-negative
  const TreeNodeIndex k = node.size();
  node.push_back (TreeNode());
  node[k].parent = -1;
  node[k].child.push_back (i);
  node[k].child.push_back (j);
  node[i].parent = k;
  node[i].d = max (0., d/2);
  node[j].parent = k;
  node[j].d = max (0., d/2);
  
  const string s = toString();
  LogThisAt(5,"Neighbor-joining tree: " << s << endl);
  parse (s);  // to ensure consistency (i.e. serializing & deserializing will not change node indices)
}

void Tree::buildByNeighborJoining (const vguard<FastSeq>& seq, const vguard<vguard<TreeBranchLength> >& distanceMatrix) {
  vguard<string> nodeName;
  nodeName.reserve (seq.size());
  for (const auto& s : seq)
    nodeName.push_back (s.name);
  return buildByNeighborJoining (nodeName, distanceMatrix);
}

void Tree::buildByUPGMA (const vguard<string>& nodeName, const vguard<vguard<TreeBranchLength> >& distanceMatrix) {
  // check that there are more than 2 nodes
  Assert (nodeName.size() >= 2, "Fewer than 2 nodes; can't make a binary tree");
  // clear the existing tree
  node.clear();
  // copy distance matrix
  vguard<vguard<TreeBranchLength> > dist = distanceMatrix;
  // estimate tree by UPGMA
  // first, initialise the list of active nodes
  set<TreeNodeIndex> activeNodes;
  for (TreeNodeIndex n = 0; n < (int) nodeName.size(); ++n) {
    activeNodes.insert (n);
    node.push_back (TreeNode());
    node.back().name = nodeName[n];
    node.back().parent = -1;
  }
  // main loop
  vguard<TreeBranchLength> nodeHeight (nodeName.size(), 0);
  while (true)
    {
      // get number of active nodes
      const int nActiveNodes = activeNodes.size();
      // loop exit test
      if (nActiveNodes == 2) break;
      Assert (nActiveNodes > 2, "Fewer than 2 nodes left -- should never get here");
      // find closest two nodes
      bool isFirstPair = true;
      TreeBranchLength minDist = 0;
      TreeNodeIndex min_i = -1, min_j = -1;
      set<TreeNodeIndex>::const_iterator node_i_p, node_j_p, activeNodes_end_p;
      activeNodes_end_p = activeNodes.end();
      for (node_i_p = activeNodes.begin(); node_i_p != activeNodes_end_p; ++node_i_p)
	for (node_j_p = node_i_p, ++node_j_p; node_j_p != activeNodes_end_p; ++node_j_p) {
	  const TreeBranchLength d = dist [*node_i_p] [*node_j_p];
	  if (isFirstPair || d < minDist)
	    {
	      min_i = *node_i_p;
	      min_j = *node_j_p;
	      minDist = d;
	      isFirstPair = false;
	    }
	}
      // nodes min_i and min_j are neighbors -- join them with new index k
      // first, calculate new distances
      const TreeNodeIndex k = nodes();
      dist.push_back (vguard<TreeBranchLength> (k + 1));
      dist[k][k] = 0;
      const TreeBranchLength d_ij = dist[min_i][min_j];
      nodeHeight.push_back (max (nodeHeight[min_i] + minBranchLength,
				 max (nodeHeight[min_j] + minBranchLength,
				      (nodeHeight[min_i] + nodeHeight[min_j] + d_ij) / 2)));
      const TreeBranchLength d_ik = nodeHeight[k] - nodeHeight[min_i];
      const TreeBranchLength d_jk = nodeHeight[k] - nodeHeight[min_j];
      for (TreeNodeIndex m = 0; m < k; ++m)
	dist[m].push_back (dist[k][m] = (dist[min_i][m] + dist[min_j][m]) / 2);
      dist[min_i][k] = dist[k][min_i] = d_ik;
      dist[min_j][k] = dist[k][min_j] = d_jk;
      // now update the Tree
      node.push_back (TreeNode());
      node[k].child.push_back (min_i);
      node[k].child.push_back (min_j);
      node[min_i].parent = k;
      node[min_i].d = max (0., d_ik);
      node[min_j].parent = k;
      node[min_j].d = max (0., d_jk);
      LogThisAt(7,"Joining nodes " << min_i << " and " << min_j << " to common ancestor " << k << " (branch lengths: " << k << "->" << min_i << " = " << d_ik << ", " << k << "->" << min_j << " = " << d_jk << ")" << endl);
      activeNodes.erase (min_i);
      activeNodes.erase (min_j);
      activeNodes.insert (k);
    }
  // make the root node
  set<TreeNodeIndex>::iterator iter = activeNodes.begin();
  const TreeNodeIndex i = *iter;
  const TreeNodeIndex j = *++iter;
  const TreeNodeIndex k = node.size();
  nodeHeight.push_back (max (nodeHeight[i] + minBranchLength,
			     max (nodeHeight[j] + minBranchLength,
				  (nodeHeight[i] + nodeHeight[j] + dist[i][j]) / 2)));
  node.push_back (TreeNode());
  node[k].parent = -1;
  node[k].child.push_back (i);
  node[k].child.push_back (j);
  node[i].parent = k;
  node[i].d = max (0., nodeHeight[k] - nodeHeight[i]);
  node[j].parent = k;
  node[j].d = max (0., nodeHeight[k] - nodeHeight[j]);

  const string s = toString();
  LogThisAt(5,"UPGMA tree: " << s << endl);
  parse (s);  // to ensure consistency (i.e. serializing & deserializing will not change node indices)

  assertUltrametric();
}

void Tree::buildByUPGMA (const vguard<FastSeq>& seq, const vguard<vguard<TreeBranchLength> >& distanceMatrix) {
  vguard<string> nodeName;
  nodeName.reserve (seq.size());
  for (const auto& s : seq)
    nodeName.push_back (s.name);
  return buildByUPGMA (nodeName, distanceMatrix);
}

string Tree::seqName (TreeNodeIndex n) const {
  string s = nodeName(n);
  if (s.size() == 0) {
    vguard<string> cs;
    for (auto c : node[n].child) {
      ostringstream o;
      o.unsetf(std::ios_base::floatfield);
      o << seqName(c) << ':' << branchLength(c);
      cs.push_back (o.str());
    }
    s = "(" + join(cs,",") + ")";
  }
  return s;
}

string Tree::pairParentName (const string& lChildName, double lTime, const string& rChildName, double rTime) {
  ostringstream o;
  o.unsetf(std::ios_base::floatfield);
  o << "(" << lChildName << ":" << lTime << "," << rChildName << ":" << rTime << ")";
  return o.str();
}

bool Tree::allNodesNamed() const {
  for (auto& n : node)
    if (n.name.empty())
      return false;
  return true;
}

void Tree::assertAllNodesNamed() const {
  Assert (allNodesNamed(), "Tree has unnamed nodes");
}
  
bool Tree::nodesMatchSeqs (const vguard<FastSeq>& seq) const {
  if (seq.size() != nodes())
    return false;
  set<string> names;
  for (size_t n = 0; n < seq.size(); ++n) {
    if (seq[n].name != node[n].name)
      return false;
    if (names.count (seq[n].name))
      return false;
    names.insert (seq[n].name);
  }
  return true;
}

void Tree::assertNodesMatchSeqs (const vguard<FastSeq>& seq) const {
  Assert (seq.size() == nodes(), "Number of sequences doesn't match number of nodes in tree");
  set<string> names;
  for (size_t n = 0; n < seq.size(); ++n) {
    Assert (seq[n].name == node[n].name, "Name of sequence #%d (%s) does not match node #%d (%s)", n, seq[n].name.c_str(), n, node[n].name.c_str());
    Assert (names.count (seq[n].name) == 0, "Duplicate sequence name %s", seq[n].name.c_str());
    names.insert (seq[n].name);
  }
  Assert (nodesMatchSeqs(seq), "Node names and sequence names don't match");  // should be redundant after previous tests...
}

bool Tree::seqNamesBijective (const vguard<FastSeq>& seq) const {
  if (!allNodesNamed())
    return false;
  if (seq.size() != nodes())
    return false;
  map<string,size_t> name2seq;
  for (size_t n = 0; n < seq.size(); ++n) {
    if (name2seq.find (seq[n].name) != name2seq.end())
      return false;
    name2seq[seq[n].name] = n;
  }
  for (TreeNodeIndex n = 0; n < nodes(); ++n)
    if (name2seq.find (seqName(n)) == name2seq.end())
      return false;
  return true;
}

void Tree::reorderSeqs (vguard<FastSeq>& seq) const {
  Assert (seq.size() == nodes(), "Number of sequences doesn't match number of nodes in tree");
  map<string,size_t> name2seq;
  for (size_t n = 0; n < seq.size(); ++n) {
    Assert (name2seq.find (seq[n].name) == name2seq.end(), "Duplicate sequence name: %s", seq[n].name.c_str());
    name2seq[seq[n].name] = n;
  }
  vguard<size_t> new2old (seq.size()), old2new (seq.size());
  for (TreeNodeIndex n = 0; n < nodes(); ++n) {
    Assert (name2seq.find (seqName(n)) != name2seq.end(), "Tree node %s is absent from sequence dataset", seqName(n).c_str());
    const size_t old_n = name2seq[seqName(n)];
    new2old[n] = old_n;
    old2new[old_n] = n;
  }

  for (size_t n = 0; n < new2old.size(); ++n) {
    const size_t o = new2old[n], m = old2new[n];
    swap (seq[n], seq[o]);
    swap (old2new[n], old2new[o]);
    swap (new2old[n], new2old[m]);
  }
}

void Tree::assignInternalNodeNames (const char* prefix) {
  set<string> names;
  for (const auto& n : node)
    if (n.name.size()) {
      Assert (names.count(n.name) == 0, "Duplicate tree node name: %s", n.name.c_str());
      names.insert (n.name);
    }
  for (size_t i = 0; i < node.size(); ++i)
    if (node[i].name.empty()) {
      string nn = string(prefix) + to_string(i+1);
      while (names.count(nn))
	nn = string("_") + nn;
      node[i].name = nn;
    }
}

void Tree::assignInternalNodeNames (vguard<FastSeq>& seq, const char* prefix) {
  reorderSeqs (seq);  // make sure that nodes match rows
  assignInternalNodeNames (prefix);
  for (size_t n = 0; n < nodes(); ++n)
    seq[n].name = seqName(n);
}

set<TreeNodeIndex> Tree::nodeAndAncestors (TreeNodeIndex node) const {
  set<TreeNodeIndex> a;
  while (node >= 0) {
    a.insert (node);
    node = parentNode(node);
  } 
  return a;
}

set<TreeNodeIndex> Tree::nodeAndDescendants (TreeNodeIndex node) const {
  const auto d = rerootedPreorderSort (node, parentNode(node));
  return set<TreeNodeIndex> (d.begin(), d.end());
}

TreeNodeIndex Tree::mostRecentCommonAncestor (TreeNodeIndex node1, TreeNodeIndex node2) const {
  const auto anc1 = nodeAndAncestors(node1);
  while (!anc1.count(node2) && node2 >= 0)
    node2 = parentNode(node2);
  return node2;
}

vguard<TreeNodeIndex> Tree::rerootedPreorderSort (TreeNodeIndex newRoot, TreeNodeIndex parentOfRoot) const {
  vguard<TreeNodeIndex> pre;
  std::function<void(TreeNodeIndex,TreeNodeIndex)> visit;
  visit = [&] (TreeNodeIndex node, TreeNodeIndex parent) -> void {
    pre.push_back (node);
    const auto kids = rerootedChildren (node, parent);
    for (auto kid : kids)
      visit (kid, node);
  };
  visit (newRoot, parentOfRoot);
  return pre;
}

vguard<TreeNodeIndex> Tree::rerootedParent (TreeNodeIndex newRoot) const {
  vguard<TreeNodeIndex> newParent (nodes());
  std::function<void(TreeNodeIndex,TreeNodeIndex)> visit;
  visit = [&] (TreeNodeIndex node, TreeNodeIndex parent) -> void {
    newParent[node] = parent;
    const auto kids = rerootedChildren (node, parent);
    for (auto kid : kids)
      visit (kid, node);
  };
  visit (newRoot, -1);
  return newParent;
}

vguard<TreeNodeIndex> Tree::preorderSort() const {
  vguard<TreeNodeIndex> roots, s;
  for (TreeNodeIndex n = 0; n < nodes(); ++n)
    if (parentNode(n) < 0)
      roots.push_back (n);
  Assert (roots.size(), "Couldn't find root");
  for (auto r: roots) {
    const auto rs = rerootedPreorderSort (r);
    s.insert (s.end(), rs.begin(), rs.end());
  }
  return s;
}

vguard<TreeNodeIndex> Tree::postorderSort() const {
  auto s = preorderSort();
  reverse (s.begin(), s.end());
  return s;
}

TreeNodeIndex Tree::closestLeaf (TreeNodeIndex node, TreeNodeIndex parent) const {
  const auto newParent = rerootedParent (parent < 0 ? node : parent);
  auto post = rerootedPreorderSort (node, parent);
  reverse (post.begin(), post.end());
  vguard<TreeNodeIndex> closest (nodes(), -1);
  vguard<double> dist (nodes());
  for (auto n : post) {
    if (isLeaf(n)) {
      closest[n] = n;
      dist[n] = 0;
    } else {
      closest[n] = -1;
      for (auto c : rerootedChildren(n,newParent[n])) {
	const TreeBranchLength d = dist[c] + branchLength(n,c);
	if (closest[n] < 0 || d < dist[n]) {
	  closest[n] = closest[c];
	  dist[n] = d;
	}
      }
    }
  }
  return closest[node];
}

Tree Tree::reorderNodes (const vguard<TreeNodeIndex>& newOrder) const {
  Tree newTree;
  newTree.node.reserve (nodes());
  vguard<TreeNodeIndex> old2new (nodes(), -1);
  for (auto oldIdx : newOrder) {
    old2new[oldIdx] = newTree.node.size();
    newTree.node.push_back (node[oldIdx]);
  }
  for (auto& n : newTree.node) {
    if (n.parent >= 0)
      n.parent = old2new[n.parent];
    for (auto& c : n.child)
      c = old2new[c];
  }
  return newTree;
}

vguard<TreeBranchLength> Tree::distanceFrom (TreeNodeIndex node) const {
  vguard<TreeBranchLength> dist (nodes());
  const vguard<TreeNodeIndex> parent = rerootedParent (node);
  for (auto n : rerootedPreorderSort(node)) {
    const auto p = parent[n];
    dist[n] = p < 0 ? 0 : max(0.,branchLength(p,n)) + dist[p];
  }
  return dist;
}

vguard<TreeBranchLength> Tree::distanceFromRoot() const {
  return distanceFrom (root());
}

void Tree::detach (TreeNodeIndex n) {
  if (node[n].parent >= 0) {
    vguard<TreeNodeIndex>& oldChild = node[node[n].parent].child;
    vguard<TreeNodeIndex> newChild;
    for (auto c : oldChild)
      if (c != n)
	newChild.push_back (c);
    oldChild.swap (newChild);
    node[n].parent = -1;
  }
}

void Tree::setParent (TreeNodeIndex n, TreeNodeIndex p, TreeBranchLength d) {
  detach (n);
  node[n].parent = p;
  node[n].d = d;
  if (p >= 0)
    node[p].child.push_back (n);
}

bool Tree::hasChildren() const {
  return nodes() > 1;
}

bool Tree::hasGrandchildren() const {
  for (TreeNodeIndex n = 0; n < root(); ++n)
    if (parentNode(n) != root())
      return true;
  return false;
}
