#ifndef PROFILE_INCLUDED
#define PROFILE_INCLUDED

#include <gsl/gsl_matrix.h>
#include "fastseq.h"
#include "alignpath.h"
#include "logsumexp.h"
#include "model.h"

typedef size_t ProfileStateIndex;
typedef size_t ProfileTransitionIndex;

struct ProfileTransition {
  ProfileStateIndex src, dest;
  LogProb lpTrans;
  AlignPath alignPath;
  IndelCounts indelCounts;
  ProfileTransition();
};

struct ProfileState {
  string name;  // for debugging only
  map<string,string> meta;  // for debugging only
  vguard<ProfileTransitionIndex> in, nullOut, absorbOut;
  vguard<LogProb> lpAbsorb;
  AlignPath alignPath;
  map<AlignRowIndex,SeqIdx> seqCoords;
  ProfileState (AlphTok alphSize = 0);
  inline bool isNull() const { return lpAbsorb.empty(); }
};

struct Profile {
  AlphTok alphSize;
  string name;  // for debugging only
  map<string,string> meta;  // for debugging only
  vguard<ProfileState> state;
  vguard<ProfileTransition> trans;
  Profile() { }
  Profile (AlphTok alphSize) : alphSize(alphSize) { }
  Profile (const string& alphabet, const FastSeq& seq, AlignRowIndex rowIndex);
  ProfileStateIndex size() const { return state.size(); }
  Profile leftMultiply (gsl_matrix* sub) const;
  const ProfileState& start() const { return state.front(); }
  const ProfileState& end() const { return state.back(); }
  const ProfileTransition* getTrans (ProfileStateIndex src, ProfileStateIndex dest) const;
  LogProb calcSumPathAbsorbProbs (const vector<LogProb>& input, const char* tag = "cumLogProb");
  void writeJson (ostream& out) const;
  string toJson() const;
};

#endif /* PROFILE_INCLUDED */
