#ifndef SAMPLER_INCLUDED
#define SAMPLER_INCLUDED

#include <iomanip>
#include "model.h"
#include "tree.h"
#include "fastseq.h"
#include "forward.h"
#include "logger.h"

struct SimpleTreePrior {
  double populationSize;
  SimpleTreePrior() : populationSize(1) { }
  double coalescenceRate (int lineages) const;
  LogProb treeLogLikelihood (const Tree& tree) const;
};

struct TreeAlignFuncs {
  typedef vguard<vguard<vguard<LogProb> > > PosWeightMatrix;  // pwm[pos][cpt][tok]

  // TreeAlignFuncs::History
  struct History {
    vguard<FastSeq> gapped;
    Tree tree;
    History reorder (const vguard<TreeNodeIndex>& newOrder) const;
    void assertNamesMatch() const;
  };

  static map<TreeNodeIndex,PosWeightMatrix> getConditionalPWMs (const RateModel& model, const Tree& tree, const vguard<FastSeq>& gapped, const map<TreeNodeIndex,TreeNodeIndex>& exclude, const set<TreeNodeIndex>& fillUpNodes, const set<TreeNodeIndex>& fillDownNodes, bool normalize = true);

  static string branchConditionalDump (const RateModel& model, const Tree& tree, const vguard<FastSeq>& gapped, TreeNodeIndex parent, TreeNodeIndex node);

  static PosWeightMatrix preMultiply (const PosWeightMatrix& child, const vguard<LogProbModel::LogProbMatrix>& submat);

  static AlignPath cladePath (const AlignPath& path, const Tree& tree, TreeNodeIndex cladeRoot, TreeNodeIndex cladeRootParent, TreeNodeIndex exclude = -1);
  static AlignPath pairPath (const AlignPath& path, TreeNodeIndex node1, TreeNodeIndex node2);
  static AlignPath triplePath (const AlignPath& path, TreeNodeIndex lChild, TreeNodeIndex rChild, TreeNodeIndex parent);
  static AlignPath branchPath (const AlignPath& path, const Tree& tree, TreeNodeIndex node);

  static bool subpathUngapped (const AlignPath& path, const vguard<TreeNodeIndex>& nodes);
  
  static set<TreeNodeIndex> allNodes (const Tree& tree);
  static set<TreeNodeIndex> allExceptNodeAndAncestors (const Tree& tree, TreeNodeIndex node);
  static set<TreeNodeIndex> nodeAndAncestors (const Tree& tree, TreeNodeIndex node);
  static set<TreeNodeIndex> nodesAndAncestors (const Tree& tree, TreeNodeIndex node1, TreeNodeIndex node2);

  static vguard<SeqIdx> getGuideSeqPos (const AlignPath& path, AlignRowIndex row, AlignRowIndex guideRow);

  static LogProb logBranchPathLikelihood (const ProbModel& probModel, const AlignPath& path, TreeNodeIndex parent, TreeNodeIndex child);
  static double rootExtProb (const RateModel& model) { return model.insExtProb; }
  static vguard<LogProb> calcInsProbs (const PosWeightMatrix& child, const vguard<LogProbModel::LogProbVector>& insvec, const vguard<LogProb>& logCptWeight);

  static LogProb rootLogLikelihood (const RateModel& model, const History& history);
  static LogProb indelLogLikelihood (const RateModel& model, const History& history);
  static LogProb substLogLikelihood (const RateModel& model, const History& history);

  static LogProb logLikelihood (const SimpleTreePrior& treePrior, const RateModel& model, const History& history, const char* suffix = "");
  static LogProb logLikelihood (const RateModel& model, const History& history, const char* suffix = "");

  // TreeAlignFuncs::SparseDPMatrix
  template <size_t CellStates>
  class SparseDPMatrix {
  protected:
    struct CellCoords {
      SeqIdx xpos, ypos;
      unsigned int state;
      CellCoords (SeqIdx xpos, SeqIdx ypos, unsigned int state)
	: xpos(xpos), ypos(ypos), state(state)
      { }
      bool operator< (const CellCoords& c) const
      { return xpos == c.xpos ? ypos == c.ypos ? state < c.state : ypos < c.ypos : xpos < c.xpos; }
      bool operator== (const CellCoords& c) const
      { return xpos == c.xpos && ypos == c.ypos && state == c.state; }
      string toString() const {
	return string("(") + to_string(xpos) + "," + to_string(ypos) + "," + to_string(state) + ")";
      }
    };

    struct XYCell {
      LogProb lp[CellStates];
      XYCell() {
	for (size_t s = 0; s < CellStates; ++s)
	  lp[s] = -numeric_limits<double>::infinity();
      }

      template<class State>
      LogProb& operator() (State s) { return lp[(unsigned int) s]; }

      template<class State>
      LogProb operator() (State s) const { return lp[(unsigned int) s]; }
    };

  private:
    const GuideAlignmentEnvelope& env;
    const vguard<SeqIdx>& xEnvPos;
    const vguard<SeqIdx>& yEnvPos;

  protected:
    const SeqIdx xSize, ySize;

  private:
    vguard<map<SeqIdx,XYCell> > cellStorage;  // partial Forward sums by cell
    XYCell emptyCell;  // always -inf
    
  public:
    LogProb lpEnd;

    // cell accessors
    inline XYCell& xyCell (SeqIdx xpos, SeqIdx ypos) { return cellStorage[xpos][ypos]; }
    inline const XYCell& xyCell (SeqIdx xpos, SeqIdx ypos) const {
      const auto& column = cellStorage[xpos];
      auto iter = column.find(ypos);
      return iter == column.end() ? emptyCell : iter->second;
    }

    template<class State>
    inline LogProb& cell (SeqIdx xpos, SeqIdx ypos, State state)
    { return cellStorage[xpos][ypos].lp[(unsigned int) state]; }

    template<class State>
    inline LogProb cell (SeqIdx xpos, SeqIdx ypos, State state) const
    {
      const auto& column = cellStorage[xpos];
      auto iter = column.find(ypos);
      return iter == column.end()
	? -numeric_limits<double>::infinity()
	: iter->second.lp[(unsigned int) state];
    }

    inline LogProb cell (const CellCoords& coords) const {
      if (coords.state == CellStates)
	return coords.xpos == xSize - 1 && coords.ypos == ySize - 1 ? lpEnd : -numeric_limits<double>::infinity();
      Assert (coords.xpos >= 0 && coords.xpos < xSize, "xpos out of range");
      Assert (coords.ypos >= 0 && coords.ypos < ySize, "ypos out of range");
      Assert (coords.state >= 0 && coords.state < CellStates, "State out of range");
      return cell (coords.xpos, coords.ypos, coords.state);
    }
    
    inline LogProb& lpStart() { return cell(0,0,0); }
    inline const LogProb lpStart() const { return cell(0,0,0); }

    inline bool inEnvelope (SeqIdx xpos, SeqIdx ypos) const {
      return xpos == 0 || ypos == 0 || xpos == xSize-1 || ypos == ySize-1
	|| env.inRange (xEnvPos[xpos], yEnvPos[ypos]);
    }
    
    // constructor
    SparseDPMatrix (const GuideAlignmentEnvelope& env, const vguard<SeqIdx>& xEnvPos, const vguard<SeqIdx>& yEnvPos)
      : env(env),
	xEnvPos(xEnvPos),
	yEnvPos(yEnvPos),
	xSize(xEnvPos.size()),
	ySize(yEnvPos.size()),
	cellStorage(xSize),
	lpEnd(-numeric_limits<double>::infinity())
    { }

    // output
    void writeToLog (int logLevel) const {
      ostringstream out;
      out << setw(5) << "x" << setw(5) << "y";
      for (size_t s = 0; s < CellStates; ++s)
	out << setw(8) << s;
      LogStream(logLevel,out.str() << endl);
      for (SeqIdx xpos = 0; xpos < xSize; ++xpos)
	for (SeqIdx ypos = 0; ypos < ySize; ++ypos)
	  if (inEnvelope (xpos, ypos)) {
	    out.str("");
	    out.clear();
	    out << setw(5) << xpos << setw(5) << ypos;
	    for (size_t s = 0; s < CellStates; ++s)
	      out << setw(10) << scientific << setprecision(2) << cell(xpos,ypos,s);
	    LogStream(logLevel,out.str() << endl);
	  }
    }
  };

  // TreeAlignFuncs::BranchMatrixBase
  class BranchMatrixBase : public SparseDPMatrix<3> {
  public:
    typedef ProbModel::State State;

    const RateModel& model;
    const ProbModel probModel;
    const LogProbModel logProbModel;

    LogProb mm, mi, md, me, im, ii, id, ie, dm, dd, de;

    AlignRowIndex xRow, yRow;
    const PosWeightMatrix& xSeq;
    const PosWeightMatrix ySub;
    const vguard<LogProb> yEmit;
   
    BranchMatrixBase (const RateModel& model, const PosWeightMatrix& xSeq, const PosWeightMatrix& ySeq, TreeBranchLength dist, const GuideAlignmentEnvelope& env, const vguard<SeqIdx>& xEnvPos, const vguard<SeqIdx>& yEnvPos, AlignRowIndex xRow, AlignRowIndex yRow);

    // helper methods
    static void getColumn (const CellCoords& coords, bool& xUngapped, bool& yUngapped);
    
    inline LogProb logMatch (SeqIdx xpos, SeqIdx ypos) const {
      return logInnerProduct (xSeq[xpos-1], ySub[ypos-1]);
    }

    LogProb lpTrans (State src, State dest) const;
    LogProb lpEmit (const CellCoords& coords) const;
    LogProb logPathProb (const AlignPath& path) const;
  };
};

struct Sampler : TreeAlignFuncs {
  typedef DPMatrix::random_engine random_engine;
  
  // Sampler::BranchMatrix
  class BranchMatrix : public BranchMatrixBase {
  public:
    BranchMatrix (const RateModel& model, const PosWeightMatrix& xSeq, const PosWeightMatrix& ySeq, TreeBranchLength dist, const GuideAlignmentEnvelope& env, const vguard<SeqIdx>& xEnvPos, const vguard<SeqIdx>& yEnvPos, AlignRowIndex xRow, AlignRowIndex yRow);

    AlignPath sample (random_engine& generator) const;
    LogProb logPostProb (const AlignPath& path) const;
  };

  // Sampler::SiblingMatrix
  struct SiblingMatrix : public SparseDPMatrix<11> {
    enum State { SSS = 0,
		 IMM = 0, IMD = 1, IDM = 2, IDD = 3,
		 WWW = 4, WWX = 5, WXW = 6,
		 IMI = 7, SSI = 7,
		 IIW = 8, SIW = 8,
		 IDI = 9, IIX = 10,
		 EEE = 11 };

    const RateModel& model;
    const ProbModel lProbModel, rProbModel;
    const LogProbModel lLogProbModel, rLogProbModel;
    const vguard<vguard<LogProb> > logRoot;  // log(cptWeight) factored in
    
    // Transition log-probabilities.
    // The null cycle idd->wxx->idd is prevented by eliminating the state wxx.
    // The outgoing transitions from wxx are folded into outgoing transitions from idd.
    // The self-transition from idd is also eliminated, and factored into outgoing transitions.
    // States {sss,ssi,siw} have same outgoing transition weights as states {imm,imi,iiw}.
    // Forward fill order: {emit states}, {www,wwx,wxw}, idd.
    // (35 transitions)
    //  To:     imm      imd      idm      idd      w**      imi      iiw      idi      iix      eee
    LogProb                                     imm_www, imm_imi, imm_iiw;
    LogProb                                     imd_wwx,                            imd_iix;
    LogProb                                     idm_wxw,                   idm_idi;
    LogProb idd_imm, idd_imd, idd_idm,                                                       idd_eee;
    LogProb www_imm, www_imd, www_idm, www_idd,                                              www_eee;
    LogProb wwx_imm, wwx_imd, wwx_idm, wwx_idd,                                              wwx_eee;
    LogProb wxw_imm, wxw_imd, wxw_idm, wxw_idd,                                              wxw_eee;
    LogProb                                     imi_www, imi_imi, imi_iiw;
    LogProb                                     iiw_www,          iiw_iiw;
    LogProb                                     idi_wxw,                   idi_idi;
    LogProb                                     iix_wwx,                            iix_iix;

    // This is 1.371* faster (48/35) but 1.375* fatter (11/8) than with w** eliminated:
    // (48 transitions)
    //        To:     imm      imd      idm      idd      imi      iiw      idi      iix      eee
    //    LogProb imm_imm, imm_imd, imm_idm, imm_idd, imm_imi, imm_iiw,                   imm_eee;
    //    LogProb imd_imm, imd_imd, imd_idm, imd_idd,                            imd_iix, imd_eee;
    //    LogProb idm_imm, idm_imd, idm_idm, idm_idd,                   idm_idi,          idm_eee;
    //    LogProb idd_imm, idd_imd, idd_idm,                                              idd_eee;
    //    LogProb imi_imm, imi_imd, imi_idm, imi_idd, imi_imi, imi_iiw,                   imi_eee;
    //    LogProb iiw_imm, iiw_imd, iiw_idm, iiw_idd,          iiw_iiw,                   iiw_eee;
    //    LogProb idi_imm, idi_imd, idi_idm, idi_idd,                   idi_idi,          idi_eee;
    //    LogProb iix_imm, iix_imd, iix_idm, iix_idd,                            iix_iix, iix_eee;
    
    AlignRowIndex lRow, rRow, pRow;
    const PosWeightMatrix lSub, rSub;
    const vguard<LogProb> lEmit, rEmit;
    
    SiblingMatrix (const RateModel& model, const PosWeightMatrix& lSeq, const PosWeightMatrix& rSeq, TreeBranchLength plDist, TreeBranchLength prDist, const GuideAlignmentEnvelope& env, const vguard<SeqIdx>& lEnvPos, const vguard<SeqIdx>& rEnvPos, AlignRowIndex lRow, AlignRowIndex rRow, AlignRowIndex pRow);

    AlignPath sample (random_engine& generator) const;
    LogProb logPostProb (const AlignPath& lrpPath) const;
    PosWeightMatrix parentSeq (const AlignPath& lrpPath) const;

    // helper methods
    static State getState (State src, bool leftUngapped, bool rightUngapped, bool parentUngapped);
    static void getColumn (const CellCoords& coords, bool& leftUngapped, bool& rightUngapped, bool& parentUngapped);
    
    LogProb lpTrans (State src, State dest) const;
    LogProb lpTransElimSelfLoopIDD (State src, State dest) const;
    LogProb lpTransElimWait (State src, State dest) const;

    inline double iddSelfLoopProb() const { return Sampler::rootExtProb(model) * lProbModel.delExt * rProbModel.delExt; }
    inline LogProb iddStay() const { return log (iddSelfLoopProb()); }
    inline LogProb iddExit() const { return log (1 / (1 - iddSelfLoopProb())); }
    
    inline LogProb rootExt() const { return log (Sampler::rootExtProb (model)); }
    inline LogProb rootNoExt() const { return log (1 - Sampler::rootExtProb (model)); }

    inline LogProb lIns() const { return log (lProbModel.ins); }
    inline LogProb lDel() const { return log (lProbModel.del); }
    inline LogProb lInsExt() const { return log (lProbModel.insExt); }
    inline LogProb lDelExt() const { return log (lProbModel.delExt); }

    inline LogProb lNoIns() const { return log (1 - lProbModel.ins); }
    inline LogProb lNoDel() const { return log (1 - lProbModel.del); }
    inline LogProb lNoInsExt() const { return log (1 - lProbModel.insExt); }
    inline LogProb lNoDelExt() const { return log (1 - lProbModel.delExt); }

    inline LogProb rIns() const { return log (rProbModel.ins); }
    inline LogProb rDel() const { return log (rProbModel.del); }
    inline LogProb rInsExt() const { return log (rProbModel.insExt); }
    inline LogProb rDelExt() const { return log (rProbModel.delExt); }

    inline LogProb rNoIns() const { return log (1 - rProbModel.ins); }
    inline LogProb rNoDel() const { return log (1 - rProbModel.del); }
    inline LogProb rNoInsExt() const { return log (1 - rProbModel.insExt); }
    inline LogProb rNoDelExt() const { return log (1 - rProbModel.delExt); }

    inline LogProb logMatch (SeqIdx xpos, SeqIdx ypos) const {
      LogProb lp = -numeric_limits<double>::infinity();
      for (int cpt = 0; cpt < model.components(); ++cpt)
	log_accum_exp (lp, logInnerProduct (logRoot[cpt], lSub[xpos-1][cpt], rSub[ypos-1][cpt]));
      return lp;
    }

    LogProb lpEmit (const CellCoords& coords) const;
  };

  // Sampler::Logger
  struct Logger {
    virtual void logHistory (const History& history) = 0;
  };
  
  // Sampler::Move
  struct Move {
    enum Type { BranchAlign = 0, NodeAlign, PruneAndRegraft, NodeHeight, Rescale, TotalMoveTypes };
    Type type;
    TreeNodeIndex node, parent, leftChild, rightChild, oldGrandparent, newGrandparent, oldSibling, newSibling;  // no single type of move uses all of these
    History oldHistory, newHistory;
    LogProb logForwardProposal, logReverseProposal, logJacobian, oldLogLikelihood, newLogLikelihood, logAcceptProb;
    bool nullified;
    string samplerName, comment;
    
    Move() { }
    Move (Type type, const History& history, LogProb logLikelihood, const string& samplerName);
    
    void initNewHistory (const Tree& tree, const vguard<FastSeq>& ungapped, const AlignPath& path);
    void initNewHistory (const Tree& tree, const vguard<FastSeq>& gapped);
    void initRatio (const Sampler& sampler);
    void nullify (const char* reason);
    bool accept (random_engine& generator) const;

    static const char* typeName (Type t);
    static size_t typeNameWidth();
  };

  struct BranchAlignMove : Move {
    BranchAlignMove (const History&, LogProb, const Sampler&, random_engine&);
  };

  struct NodeAlignMove : Move {
    NodeAlignMove (const History&, LogProb, const Sampler&, random_engine&);
  };

  struct PruneAndRegraftMove : Move {
    PruneAndRegraftMove (const History&, LogProb, const Sampler&, random_engine&);
  };

  struct NodeHeightMove : Move {
    NodeHeightMove (const History&, LogProb, const Sampler&, random_engine&);
  };

  struct RescaleMove : Move {
    RescaleMove (const History&, LogProb, const Sampler&, random_engine&);
  };

  // Sampler member variables
  const RateModel& model;
  const SimpleTreePrior& treePrior;
  list<Logger*> loggers;
  vguard<double> moveRate, moveNanosecs;
  vguard<int> movesProposed, movesAccepted;
  bool useFixedGuide, sampleAncestralSeqs;
  const Alignment guide;
  map<string,AlignRowIndex> guideRowByName;
  int maxDistanceFromGuide;

  string name;
  History currentHistory, bestHistory;
  LogProb currentLogLikelihood, bestLogLikelihood;
  bool isUltrametric;
  
  // Sampler constructor
  Sampler (const RateModel& model, const SimpleTreePrior& treePrior, const vguard<FastSeq>& gappedGuide);

  // Sampler setup methods
  void addLogger (Logger& logger);
  void initialize (const History& initialHistory, const string& name);

  // Sampler sampling methods
  inline LogProb logLikelihood (const History& history, const char* prefix = "") const {
    return TreeAlignFuncs::logLikelihood (treePrior, model, history, prefix);
  }
  
  Move proposeMove (const History& oldHistory, LogProb oldLogLikelihood, random_engine& generator) const;

  void sample (random_engine& generator);
  
  static void run (vguard<Sampler>& samplers, random_engine& generator, unsigned int nSamples = 1);

  // Sampler summary methods
  string moveStats() const;
  
  // Sampler helpers
  static TreeNodeIndex randomInternalNode (const Tree& tree, random_engine& generator);
  static TreeNodeIndex randomChildNode (const Tree& tree, random_engine& generator);
  static TreeNodeIndex randomGrandchildNode (const Tree& tree, random_engine& generator);

  static vguard<TreeNodeIndex> contemporaneousNodes (const Tree& tree, const vguard<TreeBranchLength>& distanceFromRoot, TreeNodeIndex node);
  static vguard<double> nodeListWeights (size_t n);  // returns normalized distribution over [0..n-1]
  
  AlignRowIndex guideRow (const Tree& tree, TreeNodeIndex node) const;
  GuideAlignmentEnvelope makeGuide (const Tree& tree, TreeNodeIndex leaf1, TreeNodeIndex leaf2, const AlignPath& path, TreeNodeIndex node1, TreeNodeIndex node2) const;

  vguard<SeqIdx> guideSeqPos (const AlignPath& path, AlignRowIndex row, AlignRowIndex fixedGuideRow) const;
  vguard<SeqIdx> guideSeqPos (const AlignPath& path, AlignRowIndex row, AlignRowIndex variableGuideRow, AlignRowIndex fixedGuideRow) const;

  inline map<TreeNodeIndex,PosWeightMatrix> getConditionalPWMs (const Tree& tree, const vguard<FastSeq>& gapped, const map<TreeNodeIndex,TreeNodeIndex>& exclude, const set<TreeNodeIndex>& fillUpNodes, const set<TreeNodeIndex>& fillDownNodes) const {
    return TreeAlignFuncs::getConditionalPWMs (model, tree, gapped, exclude, fillUpNodes, fillDownNodes);
  }

  string sampleSeq (const PosWeightMatrix& profile, random_engine& generator) const;
  LogProb logSeqPostProb (const string& seq, const PosWeightMatrix& profile) const;

  static string profileToString (const PosWeightMatrix& profile);
};

#endif /* SAMPLER_INCLUDED */
