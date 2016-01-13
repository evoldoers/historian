#ifndef PROFILE_INCLUDED
#define PROFILE_INCLUDED

#include "fastseq.h"

typedef size_t ProfileStateIndex;
typedef size_t ProfileTransitionIndex;
typedef double LogProb;

struct ProfileTransition {
  ProfileStateIndex src, dest;
  vector<LogProb> absorb;
  ProfileTransition (AlphTok alphSize);
};

struct ProfileState {
  vector<ProfileTransitionIndex> in, out;
};

struct Profile {
  AlphTok alphSize;
  vector<ProfileState> state;
  vector<ProfileTransition> trans;

  Profile (AlphTok alphSize, const vguard<AlphTok>& seq);
};


#endif /* PROFILE_INCLUDED */
