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
  string treeFilename, seqsFilename, modelFilename, guideFilename, reconFilename;
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
  vguard<FastSeq> seqs, gappedGuide, gappedRecon;

  map<string,size_t> seqIndex;
  map<TreeNodeIndex,size_t> nodeToSeqIndex;
  vguard<string> rowName;

  AlignPath guide;
  vguard<TreeNodeIndex> closestLeaf;
  vguard<double> closestLeafDistance;

  Alignment reconstruction;
  EventCounts counts;
  
  Reconstructor();

  bool parseReconArgs (deque<string>& argvec);
  bool parseCountArgs (deque<string>& argvec);

  bool parseTreeArgs (deque<string>& argvec);
  bool parseModelArgs (deque<string>& argvec);

  void loadReconFiles();
  void loadCountFiles();

  void reconstruct();
  void count();

  void writeRecon (ostream& out) const;
  void writeCounts (ostream& out) const;
  
private:
  void loadModel();
  void loadTree();
  void buildTree();

  void seedGenerator();

  void buildReconIndices();
};

#endif /* PROGALIGN_INCLUDED */
