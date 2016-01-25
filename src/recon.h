#ifndef RECON_INCLUDED
#define RECON_INCLUDED

#include <deque>
#include "knhx.h"
#include "alignpath.h"
#include "model.h"

struct Reconstructor {
  string treeFilename, seqsFilename, modelFilename;
  size_t profileSamples, profileNodeLimit;

  RateModel model;
  ktree_t* tree;
  vguard<FastSeq> seqs;

  map<string,size_t> seqIndex;
  map<int,size_t> nodeToSeqIndex;
  map<int,vector<int> > children;
  vguard<string> rowName;
  
  Reconstructor();
  ~Reconstructor();

  bool parseReconArgs (deque<string>& argvec);
  
  Alignment reconstruct();
  Alignment loadFilesAndReconstruct();
  void buildIndices();

  string nodeName (int node) const;
  double branchLength (int node) const;
  int nodes() const;
  int parentNode (int node) const;
  bool isLeaf (int node) const;
};


#endif /* PROGALIGN_INCLUDED */
