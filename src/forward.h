#ifndef FORWARD_INCLUDED
#define FORWARD_INCLUDED

#include <random>
#include <list>
#include "pairhmm.h"
#include "profile.h"

class ForwardMatrix {
private:
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
  vguard<LogProb> absorbScratch;  // scratch space for computing absorb profiles

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
  struct EffectiveTransition {
    LogProb lpPath, lpBestAlignPath;
    AlignPath bestAlignPath;
    EffectiveTransition();
  };

  typedef list<CellCoords> Path;
  typedef mt19937 random_engine;

  const Profile& x, y;
  const Profile subx, suby;
  const PairHMM& hmm;
  const AlignRowIndex parentRowIndex;
  const AlphTok alphSize;
  const ProfileStateIndex xSize, ySize;
  const CellCoords startCell, endCell;
  LogProb lpEnd;
  const GuideAlignmentEnvelope envelope;
  vguard<int> xClosestLeafPos, yClosestLeafPos;
  int maxDistance;

  ForwardMatrix (const Profile& x, const Profile& y, const PairHMM& hmm, AlignRowIndex parentRowIndex, const GuideAlignmentEnvelope& env);

  // cell accessors
  inline XYCell& xyCell (ProfileStateIndex xpos, ProfileStateIndex ypos) { return cellStorage[xpos][ypos]; }
  inline const XYCell& xyCell (ProfileStateIndex xpos, ProfileStateIndex ypos) const {
    const auto& column = cellStorage[xpos];
    auto iter = column.find(ypos);
    return iter == column.end() ? emptyCell : iter->second;
  }

  inline LogProb& cell (ProfileStateIndex xpos, ProfileStateIndex ypos, PairHMM::State state)
  { return cellStorage[xpos][ypos].lp[state]; }
  inline const LogProb cell (ProfileStateIndex xpos, ProfileStateIndex ypos, PairHMM::State state) const
  {
    const auto& column = cellStorage[xpos];
    auto iter = column.find(ypos);
    return iter == column.end() ? -numeric_limits<double>::infinity() : iter->second.lp[state];
  }

  // traceback
  Path sampleTrace (random_engine& generator);
  Path bestTrace();  // not quite Viterbi (takes max's rather than sampling through the Forward matrix)
  AlignPath bestAlignPath();

  // profile construction
  enum EliminationStrategy { KeepAll, KeepHubsAndAbsorbers, KeepAbsorbers };
  Profile makeProfile (const set<CellCoords>& cells, EliminationStrategy strategy = KeepHubsAndAbsorbers);
  Profile sampleProfile (random_engine& generator, size_t profileSamples, size_t maxCells = 0, EliminationStrategy strategy = KeepHubsAndAbsorbers, bool includeBestTraceInProfile = true);  // maxCells=0 to unlimit
  Profile bestProfile (EliminationStrategy strategy = KeepHubsAndAbsorbers);

  // helpers
public:
  static inline LogProb logInnerProduct (const vguard<LogProb>::const_iterator& begin1,
					 const vguard<LogProb>::const_iterator& begin2,
					 const vguard<LogProb>::const_iterator& end2) {
    LogProb lip = -numeric_limits<double>::infinity();
    for (vguard<LogProb>::const_iterator iter1 = begin1, iter2 = begin2; iter2 != end2; ++iter1, ++iter2)
      lip = log_sum_exp (lip, *iter1 + *iter2);
    return lip;
  }

  static inline LogProb logInnerProduct (const vguard<LogProb>& v1, const vguard<LogProb>& v2) {
    return logInnerProduct (v1.begin(), v2.begin(), v2.end());
  }

  string cellName (const CellCoords& cell) const;
  static string ancestorName (const string& lChildName, double lTime, const string& rChildName, double rTime);

  static random_engine newRNG();
  
private:
  inline void initAbsorbScratch (ProfileStateIndex xpos, ProfileStateIndex ypos) {
    const vguard<LogProb>::const_iterator xbegin = subx.state[xpos].lpAbsorb.begin();
    const vguard<LogProb>::const_iterator ybegin = suby.state[ypos].lpAbsorb.begin();
    for (size_t n = 0; n < hmm.alphabetSize(); ++n)
      absorbScratch[n] = xbegin[n] + ybegin[n];
  }

  inline LogProb computeLogProbAbsorb (ProfileStateIndex xpos, ProfileStateIndex ypos) {
    initAbsorbScratch (xpos, ypos);
    return logInnerProduct (hmm.logRoot, absorbScratch);
  }

  bool isAbsorbing (const CellCoords& c) const;
  
  map<CellCoords,LogProb> sourceCells (const CellCoords& destCell);
  map<CellCoords,LogProb> sourceTransitions (const CellCoords& destCell);
  LogProb eliminatedLogProbInsert (const CellCoords& cell) const;

  AlignPath cellAlignPath (const CellCoords& cell) const;
  AlignPath transitionAlignPath (const CellCoords& src, const CellCoords& dest) const;
  AlignPath traceAlignPath (const Path& path) const;

  map<AlignRowIndex,SeqIdx> cellSeqCoords (const CellCoords& cell) const;

  CellCoords sampleCell (const map<CellCoords,LogProb>& cellLogProb, random_engine& generator) const;
  static CellCoords bestCell (const map<CellCoords,LogProb>& cellLogProb);
};

#endif /* FORWARD_INCLUDED */
