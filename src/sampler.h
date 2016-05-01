#ifndef SAMPLER_INCLUDED
#define SAMPLER_INCLUDED

#include "model.h"
#include "tree.h"
#include "fastseq.h"
#include "forward.h"

struct SimpleTreePrior {
  double coalescenceRate;
  LogProb treeLogLikelihood (const Tree& tree) const;
};

struct Sampler {

  typedef DPMatrix::random_engine random_engine;

  // Sampler::SparseDPMatrix
  template <unsigned int CellStates>
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

  public:
    const TokSeq& xSeq;
    const TokSeq& ySeq;
    const GuideAlignmentEnvelope& env;
    const vguard<SeqIdx>& xEnvPos;
    const vguard<SeqIdx>& yEnvPos;

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

    // constructor
    SparseDPMatrix (const TokSeq& xSeq, const TokSeq& ySeq, const GuideAlignmentEnvelope& env, const vguard<SeqIdx>& xEnvPos, const vguard<SeqIdx>& yEnvPos)
      : xSeq(xSeq), ySeq(ySeq), env(env), xEnvPos(xEnvPos), yEnvPos(yEnvPos), lpEnd(-numeric_limits<double>::infinity())
    { }
  };
  
  // Sampler::BranchMatrix
  class BranchMatrix : public SparseDPMatrix<3> {
  public:
    enum State { Start = 0,
		 Match = 0, Insert = 1, Delete = 2,
		 End = 3,
		 SourceStates = 3, DestStates = 4 };

    const RateModel& model;
    TreeBranchLength dist;

    LogProb lpTrans[SourceStates][DestStates];
    vguard<vguard<LogProb> > submat;  // log odds-ratio

    // cell accessors
    static inline AlignRowIndex xRow() { return 0; }
    static inline AlignRowIndex yRow() { return 1; }

    BranchMatrix (const RateModel& model, const TokSeq& xSeq, const TokSeq& ySeq, TreeBranchLength dist, const GuideAlignmentEnvelope& env, const vguard<SeqIdx>& xEnvPos, const vguard<SeqIdx>& yEnvPos);

    void sample (AlignPath& path, random_engine& generator) const;
    LogProb logPostProb (const AlignPath& path) const;
  };

  // Sampler::SiblingMatrix
  struct SiblingMatrix : public SparseDPMatrix<7> {
    enum State { SSS = 0, SSI = 1, SIW = 2,
		 IMM = 0, IMI = 1, IIW = 2, IDI = 3, IIX = 4,
		 IMD = 5, IDM = 6, IDD = 7,
		 EEE = 8,
		 SourceStates = 8, DestStates = 9 };

    const RateModel& model;
    const ProbModel lProbModel, rProbModel;

    LogProb lpTrans[SourceStates][DestStates];
    LogProb lpElim[SourceStates][DestStates];  // effective transitions, with IDD eliminated
    vguard<vguard<LogProb> > lrMat;
    vguard<LogProb> lIns, rIns, lInsMat, rInsMat;

    static inline AlignRowIndex lRow() { return 0; }
    static inline AlignRowIndex rRow() { return 1; }
    static inline AlignRowIndex pRow() { return 2; }

    SiblingMatrix (const RateModel& model, const TokSeq& lSeq, const TokSeq& rSeq, TreeBranchLength plDist, TreeBranchLength prDist, const GuideAlignmentEnvelope& env, const vguard<SeqIdx>& xEnvPos, const vguard<SeqIdx>& yEnvPos);

    void sample (TokSeq& pSeq, AlignPath& plrPath, random_engine& generator) const;
    LogProb logPostProb (const TokSeq& pSeq, const AlignPath& plrPath) const;
  };

  // Sampler::History
  struct History {
    vguard<FastSeq> gapped;
    Tree tree;
  };

  // Sampler::Log
  struct Log {
    virtual void log (const History& history) = 0;
  };
  
  // Sampler::Move
  struct Move {
    enum Type { SampleBranch, SampleNode, PruneAndRegraft, SampleNodeHeight, SampleAncestralResidues };
    Type type;
    TreeNodeIndex node, parent, grandparent, leftChild, rightChild, oldSibling, newSibling;  // no single type of move uses all of these
    History oldHistory, newHistory;
    LogProb logProposal, logInverseProposal, logOldLikelihood, logNewLikelihood, logHastingsRatio;

    Move (Type type, const History& history);
    bool accept (random_engine& generator) const;
  };

  struct SampleBranchMove : Move {
    SampleBranchMove (const History&, Sampler&, random_engine&);
  };

  struct SampleNodeMove : Move {
    SampleNodeMove (const History&, Sampler&, random_engine&);
  };

  struct PruneAndRegraftMove : Move {
    PruneAndRegraftMove (const History&, Sampler&, random_engine&);
  };

  struct SampleNodeHeightMove : Move {
    SampleNodeHeightMove (const History&, Sampler&, random_engine&);
  };

  struct SampleAncestralResiduesMove : Move {
    SampleAncestralResiduesMove (const History&, Sampler&, random_engine&);
  };

  // Sampler member variables
  RateModel model;
  SimpleTreePrior treePrior;
  list<Log*> logs;
  map<Move::Type,double> moveRate;
  Alignment guide;
  int maxDistanceFromGuide;
  
  // Sampler constructor
  Sampler (const RateModel& model, const SimpleTreePrior& treePrior, const vguard<FastSeq>& gappedGuide);
  
  // Sampler methods
  void addLog (Log& log);
  Move proposeMove (const History& oldState, random_engine& generator) const;
  void run (History& state, random_engine& generator, int nSamples = 1);

  // Sampler helpers
  static TreeNodeIndex randomInternalNode (const Tree& tree, random_engine& generator);
  static vguard<SeqIdx> guideSeqPos (const AlignPath& path, AlignRowIndex row, AlignRowIndex guideRow);
  TokSeq removeGapsAndTokenize (const FastSeq& gapped) const;
  static AlignPath cladePath (const AlignPath& path, const Tree& tree, TreeNodeIndex cladeRoot, TreeNodeIndex cladeRootParent);
};

#endif /* SAMPLER_INCLUDED */
