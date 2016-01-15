#ifndef PROFILE_INCLUDED
#define PROFILE_INCLUDED

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
  vector<ProfileState> state;
  vector<ProfileTransition> trans;
  Profile (AlphTok alphSize, const vguard<AlphTok>& seq, AlignRowIndex rowIndex);
  ProfileStateIndex size() const { return state.size(); }
  Profile leftMultiply (gsl_matrix sub) const;
};

#endif /* PROFILE_INCLUDED */
