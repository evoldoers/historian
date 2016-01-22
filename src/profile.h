#ifndef PROFILE_INCLUDED
#define PROFILE_INCLUDED

#include <gsl/gsl_matrix.h>
#include "fastseq.h"
#include "alignpath.h"
#include "logsumexp.h"

typedef size_t ProfileStateIndex;
typedef size_t ProfileTransitionIndex;

struct ProfileTransition {
  ProfileStateIndex src, dest;
  LogProb lpTrans;
  AlignPath alignPath;
  ProfileTransition();
};

struct ProfileState {
  vector<ProfileTransitionIndex> in, out;
  vector<LogProb> lpAbsorb;
  AlignPath alignPath;
  ProfileState (AlphTok alphSize = 0);
  inline bool isNull() const { return lpAbsorb.empty(); }
};

struct Profile {
  AlphTok alphSize;
  vector<ProfileState> state;
  vector<ProfileTransition> trans;
  Profile() { }
  Profile (AlphTok alphSize, const vguard<AlphTok>& seq, AlignRowIndex rowIndex);
  ProfileStateIndex size() const { return state.size(); }
  Profile leftMultiply (gsl_matrix* sub) const;
  const ProfileState& start() const { return state.front(); }
  const ProfileState& end() const { return state.back(); }
  const ProfileTransition* getTrans (ProfileStateIndex src, ProfileStateIndex dest) const;
};

#endif /* PROFILE_INCLUDED */
