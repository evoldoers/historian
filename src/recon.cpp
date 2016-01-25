#include <fstream>
#include "recon.h"
#include "util.h"
#include "forward.h"
#include "jsonutil.h"
#include "logger.h"

Reconstructor::Reconstructor()
  : profileSamples (100),
    profileNodeLimit (0),
    tree (NULL)
{ }

Reconstructor::~Reconstructor()
{
  if (tree)
    kn_free (tree);
}

bool Reconstructor::parseReconArgs (deque<string>& argvec) {
  if (argvec.size()) {
    const string& arg = argvec[0];
    if (arg == "-tree") {
      Require (argvec.size() > 1, "%s must have an argument", arg.c_str());
      treeFilename = argvec[1];
      argvec.pop_front();
      argvec.pop_front();
      return true;

    } else if (arg == "-model") {
      Require (argvec.size() > 1, "%s must have an argument", arg.c_str());
      modelFilename = argvec[1];
      argvec.pop_front();
      argvec.pop_front();
      return true;

    } else if (arg == "-seqs") {
      Require (argvec.size() > 1, "%s must have an argument", arg.c_str());
      seqsFilename = argvec[1];
      argvec.pop_front();
      argvec.pop_front();
      return true;
    }
  }
  
  return false;
}

Alignment Reconstructor::loadFilesAndReconstruct() {
  Require (seqsFilename.size() > 0, "Must specify sequences");
  Require (treeFilename.size() > 0, "Must specify a tree");
  Require (modelFilename.size() > 0, "Must specify an evolutionary model");

  ifstream treeFile (treeFilename);
  tree = kn_parse (JsonUtil::readStringFromStream (treeFile).c_str());

  seqs = readFastSeqs (seqsFilename.c_str());

  ifstream modelFile (modelFilename);
  ParsedJson pj (modelFile);
  model.read (pj.value);

  return reconstruct();
}

void Reconstructor::buildIndices() {
  for (size_t n = 0; n < seqs.size(); ++n) {
    Assert (seqIndex.find (seqs[n].name) == seqIndex.end(), "Duplicate sequence name %s", seqs[n].name.c_str());
    seqIndex[seqs[n].name] = n;
  }

  for (int node = 0; node < nodes(); ++node)
    if (!isLeaf(node))
      Assert (nChildren(node) == 2, "Tree is not binary: node %d has %s\nSubtree rooted at %d: %s\n", node, plural(nChildren(node),"child","children").c_str(), node, treeString(node).c_str());

  for (int node = 0; node < nodes(); ++node)
    if (isLeaf(node)) {
      Assert (nodeName(node).length() > 0, "Leaf node %d is unnamed", node);
      Assert (seqIndex.find (nodeName(node)) != seqIndex.end(), "Can't find sequence for leaf node %s", nodeName(node).c_str());
      nodeToSeqIndex[node] = seqIndex[nodeName(node)];
      rowName.push_back (nodeName(node));
    } else
      rowName.push_back (ForwardMatrix::ancestorName (rowName[getChild(node,0)], branchLength(getChild(node,0)), rowName[getChild(node,1)], branchLength(getChild(node,1))));
}

Alignment Reconstructor::reconstruct() {
  buildIndices();
  
  gsl_vector* eqm = model.getEqmProb();
  auto generator = ForwardMatrix::newRNG();

  AlignPath path;
  map<int,Profile> prof;
  for (int node = 0; node < nodes(); ++node) {
    if (isLeaf(node))
      prof[node] = Profile (model.alphabet, seqs[nodeToSeqIndex[node]], node);
    else {
      const int lChildNode = getChild(node,0);
      const int rChildNode = getChild(node,1);
      const Profile& lProf = prof[lChildNode];
      const Profile& rProf = prof[rChildNode];
      ProbModel lProbs (model, branchLength(lChildNode));
      ProbModel rProbs (model, branchLength(rChildNode));
      PairHMM hmm (lProbs, rProbs, eqm);

      LogThisAt(1,"Aligning " << lProf.name << " and " << rProf.name);

      ForwardMatrix forward (lProf, rProf, hmm, node);
      if (node == nodes() - 1) {
	path = forward.bestAlignPath();
	prof[node] = forward.bestProfile();
      } else
	prof[node] = forward.sampleProfile (generator, profileSamples, profileNodeLimit);
    }
  }

  gsl_vector_free (eqm);

  vguard<FastSeq> ungapped (nodes());
  for (int node = 0; node < nodes(); ++node) {
    if (seqIndex.find(rowName[node]) == seqIndex.end()) {
      ungapped[node].seq = string (alignPathResiduesInRow(path.at(node)), '*');
      ungapped[node].name = rowName[node];
    } else
      ungapped[node] = seqs[seqIndex.at(rowName[node])];
  }

  return Alignment (ungapped, path);
}

string Reconstructor::nodeName (int node) const {
  if (tree->node[node].name)
    return string (tree->node[node].name);
  return string();
}

double Reconstructor::branchLength (int node) const {
  return tree->node[node].d;
}

int Reconstructor::nodes() const {
  return tree->n;
}

int Reconstructor::parentNode (int node) const {
  return tree->node[node].parent;
}

bool Reconstructor::isLeaf (int node) const {
  return nChildren (node) == 0;
}

int Reconstructor::nChildren (int node) const {
  return tree->node[node].n;
}

int Reconstructor::getChild (int node, int childNum) const {
  return tree->node[node].child[childNum];
}

string Reconstructor::treeString (int root) const {
  kstring_t str;
  str.l = str.m = 0; str.s = 0;
  kn_format(tree->node, root, &str);
  string s (str.s);
  free (str.s);
  return s;
}

string Reconstructor::treeString() const {
  return treeString (tree->n - 1);
}
