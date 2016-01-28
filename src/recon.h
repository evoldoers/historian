#ifndef RECON_INCLUDED
#define RECON_INCLUDED

#include <deque>
#include "tree.h"
#include "alignpath.h"
#include "model.h"
#include "forward.h"

class Reconstructor {
public:
  string treeFilename, seqsFilename, modelFilename, guideFilename;
  size_t profileSamples, profileNodeLimit;
  int maxDistanceFromGuide;
  bool includeBestTraceInProfile;

  ForwardMatrix::random_engine generator;
  unsigned rndSeed;

  RateModel model;
  Tree tree;
  vguard<FastSeq> seqs, gapped;

  map<string,size_t> seqIndex;
  map<TreeNodeIndex,size_t> nodeToSeqIndex;
  vguard<string> rowName;

  AlignPath guide;
  vguard<TreeNodeIndex> closestLeaf;
  vguard<double> closestLeafDistance;
  
  Reconstructor();

  bool parseReconArgs (deque<string>& argvec);
  Alignment loadFilesAndReconstruct();

private:
  void buildIndices();
  Alignment reconstruct();
};


#endif /* PROGALIGN_INCLUDED */
