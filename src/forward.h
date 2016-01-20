#ifndef FORWARD_INCLUDED
#define FORWARD_INCLUDED

#include <random>
#include <list>
#include "pairhmm.h"
#include "profile.h"

class ForwardMatrix {
private:
  vector<LogProb> cellStorage;  // partial Forward sums by cell
  vector<LogProb> insx, insy;  // insert probabilities by x & y indices
  vector<LogProb> absorbScratch;  // scratch space for computing absorb profiles

public:
  struct CellCoords {
    ProfileStateIndex xpos, ypos;
    PairHMM::State state;
    CellCoords (ProfileStateIndex xpos, ProfileStateIndex ypos, PairHMM::State state)
      : xpos(xpos), ypos(ypos), state(state)
    { }
    bool operator< (const CellCoords& c) const
    { return xpos == c.xpos ? ypos == c.ypos ? state < c.state : ypos < c.ypos : xpos < c.xpos; }
  };
  typedef list<CellCoords> Path;
  typedef default_random_engine random_engine;

  const Profile& x, y;
  const Profile subx, suby;
  const PairHMM& hmm;
  const AlphTok alphSize;
  const ProfileStateIndex xSize, ySize;
  LogProb lpEnd;
  ForwardMatrix (const Profile& x, const Profile& y, const PairHMM& hmm);
  inline double& cell (ProfileStateIndex xpos, ProfileStateIndex ypos, PairHMM::State state)
  { return cellStorage[(ypos * xSize + xpos) * PairHMM::TotalStates + state]; }
  Path sampleTrace (random_engine& generator);
  Profile makeProfile (const set<CellCoords>& states);

  // helpers
private:
  static inline LogProb logInnerProduct (gsl_vector* gslv,
					 const vector<LogProb>::const_iterator& begin,
					 const vector<LogProb>::const_iterator& end) {
    LogProb lip = -numeric_limits<double>::infinity();
    size_t n = 0;
    for (vector<LogProb>::const_iterator iter = begin; iter != end; ++iter, ++n)
      lip = log_sum_exp (lip, *iter + gsl_vector_get(gslv,n));
    return lip;
  }

  static inline LogProb logInnerProduct (gsl_vector* gslv, const vector<LogProb>& stlv) {
    return logInnerProduct (gslv, stlv.begin(), stlv.end());
  }

  inline void initAbsorbScratch (ProfileStateIndex xpos, ProfileStateIndex ypos) {
    const vector<double>::const_iterator xbegin = subx.state[xpos].lpAbsorb.begin();
    const vector<double>::const_iterator ybegin = suby.state[ypos].lpAbsorb.begin();
    for (size_t n = 0; n < hmm.alphabetSize(); ++n)
      absorbScratch[n] = xbegin[n] + ybegin[n];
  }

  inline LogProb computeLogProbAbsorb (ProfileStateIndex xpos, ProfileStateIndex ypos) {
    initAbsorbScratch (xpos, ypos);
    return logInnerProduct (hmm.root, absorbScratch);
  }

  static CellCoords sampleCell (const map<CellCoords,LogProb>& cellLogProb, random_engine& generator);
};

#endif /* FORWARD_INCLUDED */
