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
  EigenCounts counts;
  AlignPath alignPath;
  ProfileTransition();
  AlignPath bestAlignPath() const;
};

struct ProfileState {
  typedef map<AlignRowIndex,SeqIdx> SeqCoords;
  string name;  // for debugging only
  map<string,string> meta;  // for debugging only
  vguard<ProfileTransitionIndex> in, nullOut, absorbOut;
  vguard<vguard<LogProb> > lpAbsorb;
  AlignPath alignPath;
  SeqCoords seqCoords;
  ProfileState();
  ProfileState (size_t components, AlphTok alphSize);
  inline bool isNull() const { return lpAbsorb.empty(); }
  inline bool isEmit() const { return !lpAbsorb.empty(); }
  inline bool isStart() const { return in.empty(); }
  inline bool isEmitOrStart() const { return isEmit() || isStart(); }
  inline bool isReady() const { return nullOut.empty(); }
  inline bool isWait() const { return absorbOut.empty(); }
  
  static void assertSeqCoordsConsistent (const SeqCoords& srcCoords, const ProfileState& dest, const AlignPath& transPath);
  static void assertSeqCoordsConsistent (const SeqCoords& srcCoords, const SeqCoords& destCoords, const AlignPath& transPath, const AlignPath& destPath);
};

struct Profile {
  AlphTok alphSize;
  size_t components;
  string name;  // for debugging only
  map<string,string> meta;  // for debugging only
  vguard<ProfileState> state;
  vguard<ProfileTransition> trans;
  map<AlignRowIndex,string> seq;
  map<ProfileStateIndex,ProfileStateIndex> equivAbsorbState;
  Profile() { }
  Profile (size_t components, AlphTok alphSize)
    : components(components), alphSize(alphSize) { }
  Profile (size_t components, const string& alphabet, const FastSeq& seq, AlignRowIndex rowIndex);
  ProfileStateIndex size() const { return state.size(); }
  Profile leftMultiply (const vguard<gsl_matrix*>& sub) const;
  const ProfileState& start() const { return state.front(); }
  const ProfileState& end() const { return state.back(); }
  const ProfileTransition* getTrans (ProfileStateIndex src, ProfileStateIndex dest) const;
  map<AlignRowIndex,char> alignColumn (ProfileStateIndex s) const;
  LogProb calcSumPathAbsorbProbs (const vguard<LogProb>& logCptWeight, const vguard<vguard<LogProb> >& logInsProb, const char* tag = "cumLogProb");
  void writeJson (ostream& out) const;
  string toJson() const;
  string tinyDescription (ProfileStateIndex s) const;  // for debugging

  void assertTransitionsConsistent() const;
  void assertSeqCoordsConsistent() const;
  void assertAllStatesWaitOrReady() const;
  void assertPathToEndExists() const;
  Profile addReadyStates() const;

  vguard<ProfileStateIndex> examplePathToEnd() const;  // proof-of-concept path to end
};

#endif /* PROFILE_INCLUDED */
