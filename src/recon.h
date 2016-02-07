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
  vguard<string> countFilenames;
  size_t profileSamples, profileNodeLimit;
  int maxDistanceFromGuide;
  bool includeBestTraceInProfile, usePosteriorsForProfile, reconstructRoot, predictAncestralSequence, accumulateCounts;
  double minPostProb;
  
  ForwardMatrix::random_engine generator;
  unsigned rndSeed;

  DiagEnvParams diagEnvParams;
  
  RateModel model;
  Tree tree;
  vguard<FastSeq> seqs, gappedGuide, gappedRecon, ancestral;

  map<string,size_t> seqIndex;
  map<TreeNodeIndex,size_t> nodeToSeqIndex;
  vguard<string> rowName;

  AlignPath guide;
  vguard<TreeNodeIndex> closestLeaf;
  vguard<double> closestLeafDistance;

  Alignment reconstruction;
  EigenCounts eigenCounts;
  EventCounts eventCounts;
  
  Reconstructor();

  bool parseReconArgs (deque<string>& argvec);
  bool parsePostArgs (deque<string>& argvec);
  bool parseCountArgs (deque<string>& argvec);
  bool parseSumArgs (deque<string>& argvec);

  bool parseTreeArgs (deque<string>& argvec);
  bool parseModelArgs (deque<string>& argvec);

  void loadSeqs();
  void loadRecon();
  void loadCounts();

  void reconstruct();
  void count();
  void fit();

  void writeRecon (ostream& out) const;
  void writeCounts (ostream& out) const;
  void writeModel (ostream& out) const;
  
private:
  void loadModel();
  void loadTree();
  void buildTree();

  void seedGenerator();

  void buildReconIndices();
};

#endif /* PROGALIGN_INCLUDED */
