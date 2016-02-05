#ifndef RECON_INCLUDED
#define RECON_INCLUDED

#include <deque>
#include "tree.h"
#include "alignpath.h"
#include "model.h"
#include "forward.h"
#include "diagenv.h"

#define DefaultProfileSamples 100
#define DefaultProfilePostProb .1
#define DefaultMaxDistanceFromGuide 10

class Reconstructor {
public:
  string treeFilename, seqsFilename, modelFilename, guideFilename;
  string treeSaveFilename, seqsSaveFilename, modelSaveFilename, guideSaveFilename;
  size_t profileSamples, profileNodeLimit;
  int maxDistanceFromGuide;
  bool includeBestTraceInProfile, usePosteriorsForProfile;
  double minPostProb;
  
  ForwardMatrix::random_engine generator;
  unsigned rndSeed;

  DiagEnvParams diagEnvParams;
  
  RateModel model;
  Tree tree;
  vguard<FastSeq> seqs, gapped;

  map<string,size_t> seqIndex;
  map<TreeNodeIndex,size_t> nodeToSeqIndex;
  vguard<string> rowName;

  AlignPath guide;
  vguard<TreeNodeIndex> closestLeaf;
  vguard<double> closestLeafDistance;

  Alignment reconstruction;
  
  Reconstructor();

  bool parseReconArgs (deque<string>& argvec);
  bool parseModelArgs (deque<string>& argvec);

  void loadReconFiles();

  void reconstruct();

private:
  void buildIndices();
};


#endif /* PROGALIGN_INCLUDED */
