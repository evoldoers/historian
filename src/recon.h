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
#define DefaultMaxEMIterations 100
#define DefaultMinEMImprovement .01

class Reconstructor {
public:
  string fastaReconFilename, treeFilename, modelFilename;
  list<string> seqFilenames, fastaGuideFilenames, nexusGuideFilenames, stockholmGuideFilenames, nexusReconFilenames, stockholmReconFilenames, countFilenames;
  string modelSaveFilename, guideSaveFilename;
  size_t profileSamples, profileNodeLimit, maxEMIterations;
  int maxDistanceFromGuide;
  bool guideAlignTryAllPairs, includeBestTraceInProfile, keepGapsOpen, usePosteriorsForProfile, reconstructRoot, predictAncestralSequence, accumulateCounts, gotPrior, useLaplacePseudocounts;
  double minPostProb, minEMImprovement;
  typedef enum { FastaFormat, NexusFormat, StockholmFormat } OutputFormat;
  OutputFormat outputFormat;
  ofstream* guideFile;
  
  ForwardMatrix::random_engine generator;
  unsigned rndSeed;

  DiagEnvParams diagEnvParams;
  
  RateModel model;

  struct Dataset {
    Tree tree;
    vguard<FastSeq> seqs, gappedGuide, gappedRecon, gappedAncestralRecon;

    map<string,size_t> seqIndex;
    map<TreeNodeIndex,size_t> nodeToSeqIndex;
    vguard<string> rowName;

    AlignPath guide;
    vguard<TreeNodeIndex> closestLeaf;
    vguard<double> closestLeafDistance;

    Alignment reconstruction;
    EigenCounts eigenCounts;

    void initGuide (const vguard<FastSeq>& gapped);
    void prepareRecon (Reconstructor& recon);
    bool hasReconstruction() const { return !gappedRecon.empty(); }
  };
  list<Dataset> datasets;
  EventCounts priorCounts, dataCounts, dataPlusPriorCounts;
  
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
  void loadSeqs (const string& seqsFilename, const string& guideFilename, const string& nexusFilename, const string& stockholmFilename);
  void loadRecon();
  void loadCounts();

  void reconstructAll();
  void countAll();

  void reconstruct (Dataset& dataset);
  void count (Dataset& dataset);
  void fit();

  void writeTreeAlignment (const Tree& tree, const vguard<FastSeq>& gapped, ostream& out, bool isReconstruction = false) const;
  void writeRecon (const Dataset& dataset, ostream& out) const;
  void writeRecon (ostream& out) const;
  void writeCounts (ostream& out) const;
  void writeModel (ostream& out) const;
  
private:
  Dataset& newDataset();
  void loadTree (Dataset& dataset);
  void buildTree (Dataset& dataset);

  Alignment makeAlignment (const Dataset& dataset, const AlignPath& path, TreeNodeIndex root) const;
  string makeAlignmentString (const Dataset& dataset, const AlignPath& path, TreeNodeIndex root, bool assignInternalNodeNames) const;
  
  void seedGenerator();
};

#endif /* PROGALIGN_INCLUDED */
