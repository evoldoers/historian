#include "profile.h"

ProfileTransition::ProfileTransition()
  : lpTrans(-numeric_limits<double>::infinity())
{ }

ProfileState::ProfileState (AlphTok alphSize)
  : lpAbsorb(alphSize,-numeric_limits<double>::infinity())
{ }

Profile::Profile (AlphTok alphSize, const vguard<AlphTok>& seq)
  : alphSize (alphSize),
    state (seq.size() + 2, ProfileState (alphSize))
{
  state.front() = state.back() = ProfileState();  // start and end are null states
  for (size_t pos = 0; pos < seq.size(); ++pos) {
    ProfileTransition trans;
    trans.src = pos;
    trans.dest = pos + 1;
    trans.lpTrans = 0;
    state[pos].out.push_back (trans);
    state[pos+1].in.push_back (trans);
  }
}
