#ifndef FORWARD_INCLUDED
#define FORWARD_INCLUDED

#include "pairhmm.h"

class ForwardMatrix {
private:
  vector<LogProb> cellStorage;
  vector<LogProb> absorbStorage;
  inline vector<double>::iterator absorbStorageBegin (ProfileStateIndex xpos, ProfileStateIndex ypos) const
  { return absorbStorage.begin() + (ypos * x.size() + xpos) * hmm.alphabetSize(); }
  inline vector<double>::iterator absorbStorageEnd (ProfileStateIndex xpos, ProfileStateIndex ypos) const
  { return absorbStorageBegin(xpos,ypos) + hmm.alphabetSize(); }
public:
  struct CellCoords {
    ProfileStateIndex xpos, ypos;
    PairHMM::State state;
  };
  typedef vector<CellCoords> Path;
  const Profile& x, y;
  const Profile subx, suby;
  const PairHMM& hmm;
  ForwardMatrix (const Profile& x, const Profile& y, const PairHMM& hmm)
    : x(x),
      y(y),
      hmm(hmm),
      cellStorage (x.size() * y.size() * PairHMM::TotalStates, -numeric_limits<double>::infinity()),
      absorbStorage (x.size() * y.size() * hmm.alphabetSize(), -numeric_limits<double>::infinity())
  { }
  inline double& cell (ProfileStateIndex xpos, ProfileStateIndex ypos, PairHMM::State state)
  { return cellStorage[(ypos * x.size() + xpos) * PairHMM::TotalStates + state]; }
  Path sampleTrace() const;
  Path bestTrace() const;  // always chooses max direction during traceback
  Profile makeProfile (const set<CellCoords>& states) const;
};


#endif /* FORWARD_INCLUDED */
