#include "recon.h"
#include "util.h"
#include "forward.h"

Reconstructor::Reconstructor()
  : profileSamples(100),
    profileNodeLimit(0)
{ }

AlignPath Reconstructor::reconstruct (ktree_t* tree, const vguard<FastSeq>& seqs, const RateModel& model) {
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
      if (node == tree->n - 1)
	path = forward.bestAlignPath();
      else
	prof[node] = forward.sampleProfile (generator, profileSamples, profileNodeLimit);
    }
  }

  gsl_vector_free (eqm);
  
  return path;
}
