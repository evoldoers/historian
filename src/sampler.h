#ifndef SAMPLER_INCLUDED
#define SAMPLER_INCLUDED

#include "model.h"
#include "tree.h"
#include "fastseq.h"

struct SimpleTreePrior {
  double coalescenceRate;
  LogProb treeLogLikelihood (const Tree& tree) const;
};

struct Sampler {

  typedef DPMatrix::random_engine random_engine;

  struct TreeAlignState {
    vguard<FastSeq> gapped;
    Tree tree;
  };

  struct Move {
    enum Type { Branch, TripleBranch, PruneRegraft, NodeHeight };
    Type type;
    TreeNodeIndex node, parent, leftChild, rightChild, oldSibling, newSibling;
    TreeAlignState oldState, newState;
    LogProb logProposal, logInverseProposal, logOldState, logNewState, logHastingsRatio;
  };

  class AlignmentMatrix {
  public:
    enum State { Start = -1, Match = 0, Insert = 1, Delete = 2, TotalStates = 3 };

  private:
    struct XYCell {
      LogProb lp[TotalStates];
      XYCell() {
	for (size_t s = 0; s < TotalStates; ++s)
	  lp[s] = -numeric_limits<double>::infinity();
      }
      LogProb& operator() (State s) { return lp[s]; }
      LogProb operator() (State s) const { return lp[s]; }
    };
    vguard<map<SeqIdx,XYCell> > cellStorage;  // partial Forward sums by cell
    XYCell emptyCell;  // always -inf

  public:
    const RateModel* model;
    FastSeq xSeq, ySeq;
    TreeBranchLength dist;

    vguard<vguard<LogProb> > submat;  // log odds-ratio
    LogProb s2m, s2i, s2d, s2e, m2m, m2i, m2d, m2e, i2m, i2i, i2d, i2e, d2m, d2i, d2d, d2e;

    const GuideAlignmentEnvelope envelope;
    LogProb lpEnd;

    // cell accessors
    inline XYCell& xyCell (SeqIdx xpos, SeqIdx ypos) { return cellStorage[xpos][ypos]; }
    inline const XYCell& xyCell (SeqIdx xpos, SeqIdx ypos) const {
      const auto& column = cellStorage[xpos];
      auto iter = column.find(ypos);
      return iter == column.end() ? emptyCell : iter->second;
    }

    inline LogProb& cell (SeqIdx xpos, SeqIdx ypos, State state)
    { return cellStorage[xpos][ypos].lp[state]; }
    inline LogProb cell (SeqIdx xpos, SeqIdx ypos, State state) const
    {
      const auto& column = cellStorage[xpos];
      auto iter = column.find(ypos);
      return iter == column.end() ? -numeric_limits<double>::infinity() : iter->second.lp[state];
    }

    inline LogProb& cell (const CellCoords& c) { return cell(c.xpos,c.ypos,c.state); }
    inline const LogProb cell (const CellCoords& c) const { return cell(c.xpos,c.ypos,c.state); }

    inline LogProb& lpStart() { return cell(0,0,Match); }
    inline const LogProb lpStart() const { return cell(0,0,Match); }

    static inline AlignRowIndex xRow() const { return 0; }
    static inline AlignRowIndex yRow() const { return 1; }

    AlignmentMatrix (const RateModel* model, const FastSeq& xSeq, const FastSeq& ySeq, TreeBranchLength dist, const GuideAlignmentEnvelope& env, const vguard<int>& xClosestLeafPos, const vguard<int>& yClosestLeafPos);

    void sample (AlignPath& path, random_engine& generator) const;
    LogProb logPostProb (const AlignPath& path) const;
  };

  struct ParentMatrix {
    enum State { SSS, SSI, SIW, IMM, IMD, IDM, IDD, IMI, IDI, IIW, IIX, EEE };

    const RateModel* model;
    FastSeq lSeq, rSeq;
    TreeBranchLength plDist, prDist;
    AlignPath lrPath;

    vguard<map<State,LogProb> > cell;

    static inline AlignRowIndex lChild() const { return 0; }
    static inline AlignRowIndex rChild() const { return 1; }
    static inline AlignRowIndex parent() const { return 2; }

    ParentMatrix (const RateModel* model, const FastSeq& lSeq, const FastSeq& rSeq, TreeBranchLength plDist, TreeBranchLength prDist, const AlignPath& lrPath);

    void sample (FastSeq& pSeq, AlignPath& plrPath, random_engine& generator) const;
    LogProb logPostProb (const FastSeq& pSeq, const AlignPath& plrPath) const;
  };
  
  RateModel model;
  SimpleTreePrior treePrior;
};

#endif /* SAMPLER_INCLUDED */
