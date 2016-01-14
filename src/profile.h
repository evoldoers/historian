#ifndef PROFILE_INCLUDED
#define PROFILE_INCLUDED

#include "fastseq.h"

typedef size_t ProfileStateIndex;
typedef double LogProb;

struct ProfileTransition {
  ProfileStateIndex src, dest;
  LogProb lpTrans;
  ProfileTransition();
};

struct ProfileState {
  vector<ProfileTransition> in, out;
  vector<LogProb> lpAbsorb;
  ProfileState (AlphTok alphSize = 0);
  inline bool isNull() const { return lpAbsorb.empty(); }
};

struct Profile {
  AlphTok alphSize;
  vector<ProfileState> state;
  Profile (AlphTok alphSize, const vguard<AlphTok>& seq);
};

#endif /* PROFILE_INCLUDED */
