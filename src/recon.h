#ifndef RECON_INCLUDED
#define RECON_INCLUDED

#include <deque>
#include "tree.h"
#include "alignpath.h"
#include "model.h"
#include "forward.h"
#include "diagenv.h"
#include "sampler.h"

#define DefaultProfileSamples 100
#define DefaultProfilePostProb .01
#define DefaultMaxDPMemoryFraction .05

#define DefaultMaxDistanceFromGuide 40
#define DefaultMaxEMIterations 100
#define DefaultMinEMImprovement .001

#define DefaultMCMCSamplesPerSeq 100

#define AncestralSequencePostProbTag "PP"

#define ReconFastAliasArgs {"-rndspan","-kmatchn","3","-band","10","-profmaxstates","1","-jc","-norefine"}

#define DefaultSimulatorRootSeqLen 100

class Reconstructor {
public:
  typedef AlignColSumProduct::ReconPostProbMap ReconPostProbMap;

  static const vguard<string> fastAliasArgs;
  
  string fastaReconFilename, treeFilename, modelFilename;
  list<string> seqFilenames, fastaGuideFilenames, nexusGuideFilenames, stockholmGuideFilenames, nexusReconFilenames, stockholmReconFilenames, countFilenames, simulatorTreeFilenames;
  string treeRoot;
  string modelSaveFilename, guideSaveFilename, dotSaveFilename, mcmcTraceFilename;
  size_t profileSamples, profileNodeLimit, maxEMIterations, mcmcSamplesPerSeq;
  int maxDistanceFromGuide, simulatorRootSeqLen;
  bool tokenizeCodons, guideAlignTryAllPairs, jukesCantorDistanceMatrix, useUPGMA, includeBestTraceInProfile, keepGapsOpen, usePosteriorsForProfile, reconstructRoot, refineReconstruction, predictAncestralSequence, reportAncestralSequenceProbability, accumulateSubstCounts, accumulateIndelCounts, gotPrior, useLaplacePseudocounts, usePosteriorsForDot, useSeparateSubPosteriorsForDot, keepDotGapsOpen, runMCMC, outputTraceMCMC, fixGuideMCMC, outputLeavesOnly;
  double minPostProb, minEMImprovement, minDotPostProb, minDotSubPostProb;
  typedef enum { FastaFormat, GappedFastaFormat, NexusFormat, StockholmFormat, NewickFormat, JsonFormat, UnknownFormat } FileFormat;
  FileFormat outputFormat;
  ofstream* guideFile;
  size_t mcmcTraceFiles;
  
  ForwardMatrix::random_engine generator;
  unsigned rndSeed;

  DiagEnvParams diagEnvParams;
  
  RateModel model;

  struct Dataset {
    string name;
    
    Tree tree;
    vguard<FastSeq> seqs, gappedGuide, gappedRecon, gappedAncestralRecon;
    ReconPostProbMap gappedAncestralReconPostProb;

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
    void clearPrep();
    bool hasReconstruction() const { return !gappedRecon.empty(); }
    bool hasAncestralReconstruction() const { return !gappedAncestralRecon.empty(); }
  };
  vguard<Dataset> datasets;
  EventCounts priorCounts, dataCounts, dataPlusPriorCounts;
  
  Reconstructor();

  bool parseModelArgs (deque<string>& argvec);
  bool parseReconArgs (deque<string>& argvec);
  bool parseSimulatorArgs (deque<string>& argvec);
  bool parseAncSeqArgs (deque<string>& argvec);
  bool parseProfileArgs (deque<string>& argvec, bool allowReconstructions);
  bool parseSamplerArgs (deque<string>& argvec);
  bool parsePremadeArgs (deque<string>& argvec);
  bool parseCountArgs (deque<string>& argvec);
  bool parseSumArgs (deque<string>& argvec);
  bool parseFitArgs (deque<string>& argvec);

  void checkUniqueSeqFile();
  void checkUniqueTreeFile();

  void setTreeFilename (const string& fn);
  void setModelFilename (const string& fn);
  
  void loadModel();
  void loadSeqs();
  void loadSeqs (const string& seqsFilename, const string& guideFilename, const string& nexusFilename, const string& stockholmFilename);
  void loadRecon();
  void loadCounts();

  void reconstructAll();
  void refineAll();
  void predictAllAncestors();
  void countAll();
  void sampleAll();

  void reconstruct (Dataset& dataset);
  void refine (Dataset& dataset);
  void predictAncestors (Dataset& dataset);
  void count (Dataset& dataset);
  void fit();

  void simulate();
  
  struct HistoryLogger : Sampler::Logger {
    Reconstructor* recon;
    ofstream* out;
    const string& name;
    HistoryLogger (Reconstructor& recon, const string& name);
    ~HistoryLogger();
    void logHistory (const Sampler::History& history);
  };

  void writeTreeAlignment (const Tree& tree, const vguard<FastSeq>& gapped, const string& name, ostream& out, bool isReconstruction = false, const ReconPostProbMap* postProb = NULL) const;
  void writeRecon (const Dataset& dataset, ostream& out) const;
  void writeRecon (ostream& out) const;
  void writeCounts (ostream& out) const;
  void writeModel (ostream& out) const;

  static FileFormat detectFormat (const string& filename);

  static int defaultMaxProfileStates();
  
private:
  Dataset& newDataset();
  void loadTree (Dataset& dataset);
  void buildTree (Dataset& dataset);

  Alignment makeAlignment (const Dataset& dataset, const AlignPath& path, TreeNodeIndex root) const;
  string makeAlignmentString (const Dataset& dataset, const AlignPath& path, TreeNodeIndex root, bool assignInternalNodeNames) const;
  
  void seedGenerator();
};

#endif /* PROGALIGN_INCLUDED */
