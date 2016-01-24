#include <fstream>
#include "recon.h"
#include "util.h"
#include "forward.h"
#include "jsonutil.h"

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

  tree = kn_parse (treeFilename.c_str());
  seqs = readFastSeqs (seqsFilename.c_str());

  ifstream modelFile (modelFilename);
  ParsedJson pj (modelFile);
  model.read (pj.value);

  return reconstruct();
}

Alignment Reconstructor::reconstruct() {
  map<string,size_t> seqIndex;
  for (size_t n = 0; n < seqs.size(); ++n) {
    Assert (seqIndex.find (seqs[n].name) == seqIndex.end(), "Duplicate sequence name %s", seqs[n].name.c_str());
    seqIndex[seqs[n].name] = n;
  }

  map<int,vector<int> > children;
  for (int node = 0; node < tree->n; ++node)
    if (tree->node[node].parent >= 0)
      children[tree->node[node].parent].push_back (node);
  for (auto& nc : children)
    Assert (nc.second.size() == 2, "Tree is not binary: node %s has %s", tree->node[nc.first].name, plural(nc.second.size(),"child","children").c_str());
  
  for (int node = 0; node < tree->n; ++node)
    if (children.find(node) == children.end())
      Assert (seqIndex.find (tree->node[node].name) != seqIndex.end(), "Can't find sequence for node %s", tree->node[node].name);

  gsl_vector* eqm = model.getEqmProb();
  auto generator = ForwardMatrix::newRNG();
  AlignPath path;
  map<int,Profile> prof;
  for (int node = 0; node < tree->n; ++node) {
    if (children.find(node) == children.end())
      prof[node] = Profile (model.alphabet, seqs[seqIndex.at (string (tree->node[node].name))], node);
    else {
      const int lChildNode = children.at(node)[0];
      const int rChildNode = children.at(node)[1];
      const FastSeq& lSeq = seqs[seqIndex.at (string (tree->node[lChildNode].name))];
      const FastSeq& rSeq = seqs[seqIndex.at (string (tree->node[rChildNode].name))];
      Profile lProf (model.alphabet, lSeq, lChildNode);
      Profile rProf (model.alphabet, rSeq, rChildNode);
      ProbModel lProbs (model, tree->node[lChildNode].d);
      ProbModel rProbs (model, tree->node[lChildNode].d);
      PairHMM hmm (lProbs, rProbs, eqm);
      ForwardMatrix forward (lProf, rProf, hmm, node);
      if (node == tree->n - 1) {
	path = forward.bestAlignPath();
	prof[node] = forward.bestProfile();
      } else
	prof[node] = forward.sampleProfile (generator, profileSamples, profileNodeLimit);
    }
  }

  gsl_vector_free (eqm);

  vguard<FastSeq> ungapped (tree->n);
  for (int node = 0; node < tree->n; ++node) {
    const string nodeName (tree->node[node].name);
    if (seqIndex.find(nodeName) == seqIndex.end()) {
      ungapped[node].seq = string (alignPathResiduesInRow(path,node), '*');
      ungapped[node].name = nodeName.size() ? nodeName : prof[node].name;
    } else
      ungapped[node] = seqs[seqIndex.at(nodeName)];
  }

  return Alignment (ungapped, path);
}
