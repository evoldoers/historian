#ifndef FORWARD_INCLUDED
#define FORWARD_INCLUDED

#include <random>
#include <queue>
#include <list>
#include "pairhmm.h"
#include "profile.h"
#include "sumprod.h"

class DPMatrix {
protected:
  struct XYCell {
    LogProb lp[PairHMM::TotalStates];
    XYCell() {
      for (size_t s = 0; s < PairHMM::TotalStates; ++s)
	lp[s] = -numeric_limits<double>::infinity();
    }
    LogProb& operator() (PairHMM::State s) { return lp[s]; }
    LogProb operator() (PairHMM::State s) const { return lp[s]; }
  };
  vguard<map<ProfileStateIndex,XYCell> > cellStorage;  // partial Forward sums by cell
  XYCell emptyCell;  // always -inf
  vguard<LogProb> insx, insy;  // insert-on-branch probabilities by x & y indices
  vguard<LogProb> rootsubx, rootsuby;  // insert-at-root-then-substitute probabilities by x & y indices
  vguard<vguard<LogProb> > absorbScratch;  // scratch space for computing absorb profiles

public:
  struct CellCoords {
    ProfileStateIndex xpos, ypos;
    PairHMM::State state;
    CellCoords() : state(PairHMM::EEE) { }
    CellCoords (ProfileStateIndex xpos, ProfileStateIndex ypos, PairHMM::State state)
      : xpos(xpos), ypos(ypos), state(state)
    { }
    bool operator< (const CellCoords& c) const
    { return xpos == c.xpos ? ypos == c.ypos ? state < c.state : ypos < c.ypos : xpos < c.xpos; }
    bool operator== (const CellCoords& c) const
    { return xpos == c.xpos && ypos == c.ypos && state == c.state; }
  };

  enum ProfilingStrategy { KeepAll = 0, CollapseChains = 1,
			   DontCountSubstEvents = 0, CountSubstEvents = 2,
			   DontCountIndelEvents = 0, CountIndelEvents = 4,
			   DontIncludeBestTrace = 0, IncludeBestTrace = 8,
			   DontKeepGapsOpen = 0, KeepGapsOpen = 16 };

  typedef list<CellCoords> Path;
  typedef mt19937 random_engine;
  static const char* random_engine_name() { return "mt19937"; }
  
  const Profile& x, y;
  const Profile subx, suby;
  const PairHMM& hmm;
  const AlphTok alphSize;
  const ProfileStateIndex xSize, ySize;
  const CellCoords startCell, endCell;
  LogProb lpEnd;
  const GuideAlignmentEnvelope envelope;
  vguard<SeqIdx> xClosestLeafPos, yClosestLeafPos;
  int maxDistance;

  DPMatrix (const Profile& x, const Profile& y, const PairHMM& hmm, const GuideAlignmentEnvelope& env);

  // cell accessors
  inline XYCell& xyCell (ProfileStateIndex xpos, ProfileStateIndex ypos) { return cellStorage[xpos][ypos]; }
  inline const XYCell& xyCell (ProfileStateIndex xpos, ProfileStateIndex ypos) const {
    const auto& column = cellStorage[xpos];
    auto iter = column.find(ypos);
    return iter == column.end() ? emptyCell : iter->second;
  }

  inline LogProb& cell (ProfileStateIndex xpos, ProfileStateIndex ypos, PairHMM::State state)
  { return cellStorage[xpos][ypos].lp[state]; }
  inline LogProb cell (ProfileStateIndex xpos, ProfileStateIndex ypos, PairHMM::State state) const
  {
    const auto& column = cellStorage[xpos];
    auto iter = column.find(ypos);
    return iter == column.end() ? -numeric_limits<double>::infinity() : iter->second.lp[state];
  }

  inline LogProb& cell (const CellCoords& c) { return cell(c.xpos,c.ypos,c.state); }
  inline const LogProb cell (const CellCoords& c) const { return cell(c.xpos,c.ypos,c.state); }

  inline LogProb& lpStart() { return cell(0,0,PairHMM::IMM); }
  inline const LogProb lpStart() const { return cell(0,0,PairHMM::IMM); }
  
  // helpers
public:
  string cellName (const CellCoords& cell) const;
  static string ancestorName (const string& lChildName, double lTime, const string& rChildName, double rTime);

  static random_engine newRNG();

  static size_t cellSize() { return sizeof(XYCell); }

  inline int components() const { return hmm.components(); }
  
protected:
  inline void initAbsorbScratch (ProfileStateIndex xpos, ProfileStateIndex ypos) {
    for (int cpt = 0; cpt < components(); ++cpt) {
      const vguard<LogProb>::const_iterator xbegin = subx.state[xpos].lpAbsorb[cpt].begin();
      const vguard<LogProb>::const_iterator ybegin = suby.state[ypos].lpAbsorb[cpt].begin();
      for (size_t n = 0; n < hmm.alphabetSize(); ++n)
	absorbScratch[cpt][n] = xbegin[n] + ybegin[n];
    }
  }

  inline LogProb computeLogProbAbsorb (ProfileStateIndex xpos, ProfileStateIndex ypos) {
    initAbsorbScratch (xpos, ypos);
    return logInnerProduct (hmm.logRoot, absorbScratch);
  }

  LogProb lpCellEmitOrAbsorb (const CellCoords& c);
  
  bool isAbsorbing (const CellCoords& c) const;
  bool changesX (const CellCoords& c) const;
  bool changesY (const CellCoords& c) const;
  
  list<CellCoords> equivAbsorbCells (const CellCoords& c) const;
  
  CellCoords sampleCell (const map<CellCoords,LogProb>& cellLogProb, random_engine& generator) const;
  static CellCoords bestCell (const map<CellCoords,LogProb>& cellLogProb);
};

class ForwardMatrix : public DPMatrix {
public:
  const AlignRowIndex parentRowIndex;
  SumProduct *sumProd;
  map<ProfileStateIndex,EigenCounts> xInsertCounts, yInsertCounts;

  struct EffectiveTransition {
    LogProb lpPath, lpBestAlignPath;
    AlignPath bestAlignPath;
    EigenCounts counts;
    EffectiveTransition();
  };
  
  ForwardMatrix (const Profile& x, const Profile& y, const PairHMM& hmm, AlignRowIndex parentRowIndex, const GuideAlignmentEnvelope& env, SumProduct* sumProd = NULL);

  // traceback
  Path sampleTrace (random_engine& generator);
  Path bestTrace();  // not quite Viterbi (takes max's rather than sampling through the Forward matrix)
  Path bestTrace (const CellCoords& end);
  AlignPath bestAlignPath();

  // profile construction
  Profile makeProfile (const set<CellCoords>& cells, ProfilingStrategy strategy = CollapseChains);
  Profile sampleProfile (random_engine& generator, size_t profileSamples, size_t maxCells = 0, ProfilingStrategy strategy = CollapseChains);  // maxCells=0 to unlimit
  Profile bestProfile (ProfilingStrategy strategy = CollapseChains);

  map<AlignRowIndex,char> getAlignmentColumn (const CellCoords& cell) const;

  void accumulateEigenCounts (EigenCounts& counts, const CellCoords& cell, SumProduct& sumProd, double weight = 1.) const;
  void accumulateCachedEigenCounts (EigenCounts& counts, const CellCoords& cell, SumProduct& sumProd, double weight = 1.);

  EigenCounts transitionEigenCounts (const CellCoords& src, const CellCoords& dest) const;
  EigenCounts cellEigenCounts (const CellCoords& cell, SumProduct& sumProd) const;
  EigenCounts cachedCellEigenCounts (const CellCoords& cell, SumProduct& sumProd);

  map<CellCoords,LogProb> sourceTransitions (const CellCoords& destCell);
  map<CellCoords,LogProb> sourceTransitionsWithoutEmitOrAbsorb (const CellCoords& destCell);

  void slowFillTest();

private:
  map<CellCoords,LogProb> sourceCells (const CellCoords& destCell);
  LogProb eliminatedLogProbInsert (const CellCoords& cell) const;

  AlignPath cellAlignPath (const CellCoords& cell) const;
  AlignPath transitionAlignPath (const CellCoords& src, const CellCoords& dest) const;
  AlignPath traceAlignPath (const Path& path) const;
  
  ProfileState::SeqCoords cellSeqCoords (const CellCoords& cell) const;
};

class BackwardMatrix : public DPMatrix {
public:
  struct CellPostProb : CellCoords {
    LogProb logPostProb;
    CellPostProb (ProfileStateIndex xpos, ProfileStateIndex ypos, PairHMM::State state, LogProb lpp)
      : CellCoords(xpos,ypos,state), logPostProb(lpp)
    { }
    bool operator< (const CellPostProb& cpp) const
    { return logPostProb < cpp.logPostProb; }
  };
  ForwardMatrix& fwd;
  
  BackwardMatrix (ForwardMatrix& fwd);

  // posterior probabilities & counts
  double cellPostProb (const CellCoords& cell) const;
  double transPostProb (const CellCoords& src, const CellCoords& dest) const;

  EigenCounts getCounts() const;

  // traceforward
  Path bestTrace (const CellCoords& start);

  // profile construction
  priority_queue<CellPostProb> cellsAbovePostProbThreshold (double minPostProb) const;
  Profile postProbProfile (double minPostProb, size_t maxCells = 0, ProfilingStrategy strategy = CollapseChains);  // maxCells=0 to unlimit
  Profile bestProfile (ProfilingStrategy strategy = CollapseChains);

  map<CellCoords,LogProb> destTransitions (const CellCoords& srcCell);

  void slowFillTest();
  void sourceDestTransTest();

private:
  map<CellCoords,LogProb> destCells (const CellCoords& srcCell);

  bool addCells (set<CellCoords>& cells, size_t maxCells, const list<CellCoords>& fwdTrace, const list<CellCoords>& backTrace, bool keepGapsOpen);
  bool addTrace (const CellCoords& cell, set<CellCoords>& cells, size_t maxCells, bool keepGapsOpen);
};

#endif /* FORWARD_INCLUDED */
