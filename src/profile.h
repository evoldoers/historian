#ifndef PROFILE_INCLUDED
#define PROFILE_INCLUDED

#include "fastseq.h"

typedef size_t ProfileStateIndex;
typedef size_t ProfileTransitionIndex;
typedef double LogProb;

struct ProfileTransition {
  ProfileStateIndex src, dest;
  LogProb lpTrans;
  vector<bool> leftAbsorbPath, rightAbsorbPath;
  ProfileTransition();
};

struct ProfileState {
  vector<ProfileTransitionIndex> in, out;
  vector<LogProb> lpAbsorb;
  bool leftAbsorb, rightAbsorb;
  ProfileState (AlphTok alphSize = 0);
  inline bool isNull() const { return lpAbsorb.empty(); }
};

struct Profile {
  vector<ProfileState> state;
  vector<ProfileTransition> trans;
  Profile (AlphTok alphSize, const vguard<AlphTok>& seq);
  ProfileStateIndex size() const { return state.size(); }
  Profile leftMultiply (gsl_matrix sub) const;
};

#endif /* PROFILE_INCLUDED */
