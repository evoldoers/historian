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
  string reconFilename, treeFilename, modelFilename;
  list<string> seqsFilenames, guideFilenames, nexusFilenames;
  string treeSaveFilename, seqsSaveFilename, modelSaveFilename, guideSaveFilename;
  vguard<string> countFilenames;
  size_t profileSamples, profileNodeLimit;
  int maxDistanceFromGuide;
  bool includeBestTraceInProfile, keepGapsOpen, usePosteriorsForProfile, reconstructRoot, predictAncestralSequence, accumulateCounts;
  double minPostProb;
  
  ForwardMatrix::random_engine generator;
  unsigned rndSeed;

  DiagEnvParams diagEnvParams;
  
  RateModel model;

  struct Dataset {
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

    void initGuide (const vguard<FastSeq>& gapped);
    void buildReconIndices();
  };
  list<Dataset> datasets;
  EventCounts eventCounts;
  
  Reconstructor();

  bool parseReconArgs (deque<string>& argvec);
  bool parsePostArgs (deque<string>& argvec);
  bool parseCountArgs (deque<string>& argvec);
  bool parseSumArgs (deque<string>& argvec);

  void checkUniqueSeqFile();
  void checkUniqueTreeFile();

  bool parseTreeArgs (deque<string>& argvec);
  bool parseModelArgs (deque<string>& argvec);

  void loadModel();
  void loadSeqs();
  void loadSeqs (const string& seqsFilename, const string& guideFilename, const string& nexusFilename);
  void loadRecon();
  void loadCounts();

  void reconstructAll();
  void countAll();

  void reconstruct (Dataset& dataset);
  void count (Dataset& dataset);
  void fit();

  void writeRecon (const Dataset& dataset, ostream& out) const;
  void writeRecon (ostream& out) const;
  void writeCounts (ostream& out) const;
  void writeModel (ostream& out) const;
  
private:
  void loadTree (Dataset& dataset);
  void buildTree (Dataset& dataset);

  void seedGenerator();
};

#endif /* PROGALIGN_INCLUDED */
