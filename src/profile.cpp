#include "profile.h"

ProfileTransition::ProfileTransition()
  : lpTrans(-numeric_limits<double>::infinity())
{ }

ProfileState::ProfileState (AlphTok alphSize)
  : lpAbsorb(alphSize,-numeric_limits<double>::infinity()),
    leftAbsorb(false),
    rightAbsorb(false)
{ }

Profile::Profile (AlphTok alphSize, const vguard<AlphTok>& seq)
  : alphSize (alphSize),
    state (seq.size() + 2, ProfileState (alphSize)),
    trans (seq.size() + 1)
{
  state.front() = state.back() = ProfileState();  // start and end are null states
  for (size_t pos = 0; pos <= seq.size(); ++pos) {
    ProfileTransition& t (trans[pos]);
    trans.src = pos;
    trans.dest = pos + 1;
    trans.lpTrans = 0;
    state[pos].out.push_back (pos);
    state[pos+1].in.push_back (pos);
  }
}
