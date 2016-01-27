#ifndef RECON_INCLUDED
#define RECON_INCLUDED

#include <deque>
#include "knhx.h"
#include "alignpath.h"
#include "model.h"

class Reconstructor {
public:
  string treeFilename, seqsFilename, modelFilename, guideFilename;
  size_t profileSamples, profileNodeLimit;
  int maxDistanceFromGuide;

  RateModel model;
  ktree_t* tree;
  vguard<FastSeq> seqs;

  map<string,size_t> seqIndex;
  map<int,size_t> nodeToSeqIndex;
  vguard<string> rowName;

  AlignPath guide;
  vguard<int> closestLeaf;
  vguard<double> closestLeafDistance;
  
  Reconstructor();
  ~Reconstructor();

  bool parseReconArgs (deque<string>& argvec);
  Alignment loadFilesAndReconstruct();

  // tree accessors
  string nodeName (int node) const;
  double branchLength (int node) const;
  int nodes() const;
  int parentNode (int node) const;
  bool isLeaf (int node) const;
  int nChildren (int node) const;
  int getChild (int node, int childNum) const;
  string treeString (int root) const;
  string treeString() const;

private:
  void buildIndices();
  Alignment reconstruct();
};


#endif /* PROGALIGN_INCLUDED */
