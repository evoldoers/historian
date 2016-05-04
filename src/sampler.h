#ifndef SAMPLER_INCLUDED
#define SAMPLER_INCLUDED

#include "model.h"
#include "tree.h"
#include "fastseq.h"
#include "forward.h"

struct SimpleTreePrior {
  double populationSize;
  double coalescenceRate (int lineages) const;
  LogProb treeLogLikelihood (const Tree& tree) const;
};

struct Sampler {

  typedef DPMatrix::random_engine random_engine;
  typedef vguard<vguard<LogProb> > PosWeightMatrix;

  // Sampler::SparseDPMatrix
  template <size_t CellStates>
  class SparseDPMatrix {
  private:
    struct XYCell {
      LogProb lp[CellStates];
      XYCell() {
	for (size_t s = 0; s < CellStates; ++s)
	  lp[s] = -numeric_limits<double>::infinity();
      }
      LogProb& operator() (unsigned int s) { return lp[s]; }
      LogProb operator() (unsigned int s) const { return lp[s]; }
    };
    vguard<map<SeqIdx,XYCell> > cellStorage;  // partial Forward sums by cell
    XYCell emptyCell;  // always -inf

    const GuideAlignmentEnvelope& env;
    const vguard<SeqIdx>& xEnvPos;
    const vguard<SeqIdx>& yEnvPos;

  public:
    LogProb lpEnd;

    // cell accessors
    inline XYCell& xyCell (SeqIdx xpos, SeqIdx ypos) { return cellStorage[xpos][ypos]; }
    inline const XYCell& xyCell (SeqIdx xpos, SeqIdx ypos) const {
      const auto& column = cellStorage[xpos];
      auto iter = column.find(ypos);
      return iter == column.end() ? emptyCell : iter->second;
    }

    inline LogProb& cell (SeqIdx xpos, SeqIdx ypos, unsigned int state)
    { return cellStorage[xpos][ypos].lp[state]; }
    inline LogProb cell (SeqIdx xpos, SeqIdx ypos, unsigned int state) const
    {
      const auto& column = cellStorage[xpos];
      auto iter = column.find(ypos);
      return iter == column.end() ? -numeric_limits<double>::infinity() : iter->second.lp[state];
    }

    inline LogProb& lpStart() { return cell(0,0,0); }
    inline const LogProb lpStart() const { return cell(0,0,0); }

    inline bool inEnvelope (SeqIdx xpos, SeqIdx ypos) const {
      return env.inRange (xEnvPos[xpos], yEnvPos[ypos]);
    }
    
    // constructor
    SparseDPMatrix (const GuideAlignmentEnvelope& env, const vguard<SeqIdx>& xEnvPos, const vguard<SeqIdx>& yEnvPos)
      : env(env), xEnvPos(xEnvPos), yEnvPos(yEnvPos), lpEnd(-numeric_limits<double>::infinity())
    { }
  };
  
  // Sampler::BranchMatrix
  class BranchMatrix : public SparseDPMatrix<3> {
  public:
    typedef ProbModel::State State;

    const RateModel& model;
    const ProbModel probModel;
    const LogProbModel logProbModel;

    LogProb mm, mi, md, me, im, ii, id, ie, dm, dd, de;
    vguard<vguard<LogProb> > submat;  // log odds-ratio

    AlignRowIndex xRow, yRow;
    const PosWeightMatrix& xSeq;
    const PosWeightMatrix ySub;
    const vguard<LogProb> yIns;
   
    BranchMatrix (const RateModel& model, const PosWeightMatrix& xSeq, const PosWeightMatrix& ySeq, TreeBranchLength dist, const GuideAlignmentEnvelope& env, const vguard<SeqIdx>& xEnvPos, const vguard<SeqIdx>& yEnvPos, AlignRowIndex xRow, AlignRowIndex yRow);

    AlignPath sample (random_engine& generator) const;
    LogProb logPostProb (const AlignPath& path) const;
  };

  // Sampler::SiblingMatrix
  struct SiblingMatrix : public SparseDPMatrix<11> {
    enum State { SSS = 0, SSI = 5, SIW = 6,
		 IMM = 0, IMD = 1, IDM = 2, IDD = 3,
		 WWW = 4, WWX = 5, WXW = 6,
		 IMI = 7, IIW = 8,
		 IDI = 9, IIX = 10,
		 EEE = 11,
		 SourceStates = 10, DestStates = 11 };

    const RateModel& model;
    const ProbModel lProbModel, rProbModel;
    const LogProbModel lLogProbModel, rLogProbModel;

    // Transition log-probabilities.
    // The null cycle idd->wxx->idd is prevented by eliminating the state wxx.
    // The outgoing paths from wxx are replaced with the following transitions:
    //  wxx->wwx->{imd,eee} is folded into idd->{imd,eee} (prevents degeneracy between wxx->wwx->www and wxx->wxw->www),
    //  wxx->wxw is replaced with idd->wxw,
    //  wxx->eee is replaced with idd->eee.
    // States {sss,ssi,siw} have same outgoing transition weights as states {imm,imi,iiw}.
    // Forward fill order: {emit states}, {wwx,wxx}, wxw, www, idd.
    // (32 transitions)
    //  To:     imm      imd      idm      idd      w**      imi      iiw      idi      iix      eee
    LogProb                                     imm_www, imm_imi, imm_iiw;
    LogProb                                     imd_wwx,                            imd_iix;
    LogProb                                     idm_wxw,                   idm_idi;
    LogProb idd_imm, idd_imd, idd_idm,          idd_wxw,                                     idd_eee;
    LogProb www_imm, www_imd, www_idm, www_idd,                                              www_eee;
    LogProb          wwx_imd,          wwx_idd, wwx_www;
    LogProb                   wxw_idm, wxw_idd, wxw_www;
    LogProb                                     imi_www, imi_imi, imi_iiw;
    LogProb                                     iiw_www,          iiw_iiw;
    LogProb                                     idi_wxw,                  idi_idi;
    LogProb                                     iix_wwx,                           iix_iix;

    // This is 1.5* faster (48/32) but 1.375* fatter (11/8) than with w** eliminated:
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
    
    vguard<vguard<LogProb> > submat;

    AlignRowIndex lRow, rRow, pRow;
    const PosWeightMatrix lSub, rSub;
    const vguard<LogProb> lIns, rIns;
    
    SiblingMatrix (const RateModel& model, const PosWeightMatrix& lSeq, const PosWeightMatrix& rSeq, TreeBranchLength plDist, TreeBranchLength prDist, const GuideAlignmentEnvelope& env, const vguard<SeqIdx>& lEnvPos, const vguard<SeqIdx>& rEnvPos, AlignRowIndex lRow, AlignRowIndex rRow, AlignRowIndex pRow);

    AlignPath sample (random_engine& generator) const;
    LogProb logPostProb (const AlignPath& plrPath) const;
    PosWeightMatrix parentSeq (const AlignPath& plrPath) const;
  };

  // Sampler::History
  struct History {
    vguard<FastSeq> gapped;
    Tree tree;
    void swapNodes (TreeNodeIndex x, TreeNodeIndex y);
  };

  // Sampler::Logger
  struct Logger {
    virtual void log (const History& history) = 0;
  };
  
  // Sampler::Move
  struct Move {
    enum Type { BranchAlign, NodeAlign, PruneAndRegraft, NodeHeight };
    Type type;
    TreeNodeIndex node, parent, leftChild, rightChild, oldGrandparent, newGrandparent, oldSibling, newSibling;  // no single type of move uses all of these
    History oldHistory, newHistory;
    LogProb logForwardProposal, logReverseProposal, oldLogLikelihood, newLogLikelihood, logHastingsRatio;

    Move (Type type, const History& history);
    void initNewHistory (const Tree& tree, const vguard<FastSeq>& ungapped, const AlignPath& path);
    void initNewHistory (const Tree& tree, const vguard<FastSeq>& gapped);
    void initRatio (const Sampler& sampler);
    bool accept (random_engine& generator) const;
  };

  struct BranchAlignMove : Move {
    BranchAlignMove (const History&, Sampler&, random_engine&);
  };

  struct NodeAlignMove : Move {
    NodeAlignMove (const History&, Sampler&, random_engine&);
  };

  struct PruneAndRegraftMove : Move {
    PruneAndRegraftMove (const History&, Sampler&, random_engine&);
  };

  struct NodeHeightMove : Move {
    NodeHeightMove (const History&, Sampler&, random_engine&);
  };

  // Sampler member variables
  RateModel model;
  SimpleTreePrior treePrior;
  list<Logger*> loggers;
  map<Move::Type,double> moveRate;
  Alignment guide;
  int maxDistanceFromGuide;
  
  // Sampler constructor
  Sampler (const RateModel& model, const SimpleTreePrior& treePrior, const vguard<FastSeq>& gappedGuide);
  
  // Sampler methods
  void addLogger (Logger& logger);
  LogProb logLikelihood (const History& history) const;
  Move proposeMove (const History& oldHistory, random_engine& generator) const;
  void run (History& state, random_engine& generator, int nSamples = 1);

  // Sampler helpers
  static TreeNodeIndex randomInternalNode (const Tree& tree, random_engine& generator);
  static TreeNodeIndex randomChildNode (const Tree& tree, random_engine& generator);
  static TreeNodeIndex randomGrandchildNode (const Tree& tree, random_engine& generator);
  static TreeNodeIndex randomContemporaneousNode (const Tree& tree, const vguard<TreeBranchLength>& distanceFromRoot, TreeNodeIndex node, random_engine& generator);

  static vguard<SeqIdx> guideSeqPos (const AlignPath& path, AlignRowIndex row, AlignRowIndex guideRow);
  map<TreeNodeIndex,PosWeightMatrix> getConditionalPWMs (const History& history, const map<TreeNodeIndex,TreeNodeIndex>& exclude) const;

  static AlignPath cladePath (const AlignPath& path, const Tree& tree, TreeNodeIndex cladeRoot, TreeNodeIndex cladeRootParent);
  static AlignPath pairPath (const AlignPath& path, TreeNodeIndex node1, TreeNodeIndex node2);
  static AlignPath triplePath (const AlignPath& path, TreeNodeIndex lChild, TreeNodeIndex rChild, TreeNodeIndex parent);
  static AlignPath branchPath (const AlignPath& path, const Tree& tree, TreeNodeIndex node);

  static LogProb logBranchPathLikelihood (const ProbModel& probModel, const AlignPath& path, TreeNodeIndex parent, TreeNodeIndex child);

  static PosWeightMatrix preMultiply (const PosWeightMatrix& child, const LogProbModel::LogProbMatrix& submat);
  static vguard<LogProb> calcInsProbs (const PosWeightMatrix& child, const LogProbModel::LogProbVector& insvec);
};

#endif /* SAMPLER_INCLUDED */
