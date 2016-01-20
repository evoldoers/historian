#ifndef FORWARD_INCLUDED
#define FORWARD_INCLUDED

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
  };
  typedef vector<CellCoords> Path;
  const Profile& x, y;
  const Profile subx, suby;
  const PairHMM& hmm;
  const size_t alphSize;
  ForwardMatrix (const Profile& x, const Profile& y, const PairHMM& hmm);
  inline double& cell (ProfileStateIndex xpos, ProfileStateIndex ypos, PairHMM::State state)
  { return cellStorage[(ypos * x.size() + xpos) * PairHMM::TotalStates + state]; }
  Path sampleTrace();
  Path bestTrace();  // always chooses max direction during traceback
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
};

#endif /* FORWARD_INCLUDED */
