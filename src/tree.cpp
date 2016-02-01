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

string Tree::nodeToString (TreeNodeIndex root) const {
  if (isLeaf(root))
    return nodeName(root);
  string s = "(";
  for (int c = 0; c < nChildren(root); ++c) {
    ostringstream d;
    d << defaultfloat << branchLength(getChild(root,c));
    if (c > 0)
      s += ",";
    s += nodeToString(getChild(root,c)) + ":" + d.str(); // to_string(branchLength(getChild(root,c)));
  }
  s += ")" + nodeName(root);
  return s;
}

string Tree::toString() const {
  return toString(root());
}

string Tree::toString (TreeNodeIndex n) const {
  return nodeToString(n) + ";";
}

string Tree::nodeName (TreeNodeIndex n) const {
  return node[n].name;
}

TreeBranchLength Tree::branchLength (TreeNodeIndex n) const {
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

TreeNodeIndex Tree::getSibling (TreeNodeIndex node) const {
  const TreeNodeIndex parent = parentNode(node);
  Assert (parent >= 0, "Attempt to find sibling of root node");
  Assert (nChildren(parent) == 2, "Attempt to find sibling in non-binary tree");
  return getChild(parent,0) == node ? getChild(parent,1) : getChild(parent,0);
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
	  LogThisAt(4,"Distance correction for node " << ni << " is " << avgDist[ni] << endl);
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
	    LogThisAt(4,"Compensated distance from node " << *node_i_p << " to node " << *node_j_p << " is " << compensatedDist << endl);
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
      LogThisAt(6,"Before Kuhner-Felsenstein:\ni=" << min_i << ", j=" << min_j << ", k=" << k << ", d_ij=" << d_ij << ", d_ik=" << d_ik << ", d_jk=" << d_jk << "\nDistances from k to other nodes: " << to_string_join(dist[k]) << endl);
      // apply Kuhner-Felsenstein correction to prevent negative branch lengths
      if (d_ik < 0)
	{
	  d_jk -= d_ik;
	  d_ik = 0;
	}
      if (d_jk < 0)
	{
	  d_ik -= d_jk;
	  d_jk = 0;
	}
      dist[min_i][k] = dist[k][min_i] = d_ik;
      dist[min_j][k] = dist[k][min_j] = d_jk;
      // now update the Tree
      node.push_back (TreeNode());
      node[k].child.push_back (min_i);
      node[k].child.push_back (min_j);
      node[min_i].parent = k;
      node[min_i].d = d_ik;
      node[min_j].parent = k;
      node[min_j].d = d_jk;
      LogThisAt(4,"Joining nodes " << min_i << " and " << min_j << " to common ancestor " << k << " (branch lengths: " << k << "->" << min_i << " = " << d_ik << ", " << k << "->" << min_j << " = " << d_jk << ")" << endl);
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
  node[k].child.push_back (i);
  node[k].child.push_back (j);
  node[i].parent = k;
  node[i].d = d/2;
  node[j].parent = k;
  node[j].d = d/2;

  const string s = toString();
  LogThisAt(4,"Neighbor-joining tree: " << s << endl);
  parse (s);  // to ensure consistency
}

void Tree::buildByNeighborJoining (const vguard<FastSeq>& seq, const vguard<vguard<TreeBranchLength> >& distanceMatrix) {
  vguard<string> nodeName;
  nodeName.reserve (seq.size());
  for (const auto& s : seq)
    nodeName.push_back (s.name);
  return buildByNeighborJoining (nodeName, distanceMatrix);
}
