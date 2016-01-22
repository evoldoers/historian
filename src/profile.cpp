#include "profile.h"

ProfileTransition::ProfileTransition()
  : lpTrans(-numeric_limits<double>::infinity())
{ }

ProfileState::ProfileState (AlphTok alphSize)
  : lpAbsorb(alphSize,-numeric_limits<double>::infinity())
{ }

Profile::Profile (AlphTok alphSize, const vguard<AlphTok>& seq, AlignRowIndex rowIndex)
  : alphSize (alphSize),
    state (seq.size() + 2, ProfileState (alphSize)),
    trans (seq.size() + 1)
{
  state.front() = state.back() = ProfileState();  // start and end are null states
  for (size_t pos = 0; pos <= seq.size(); ++pos) {
    ProfileTransition& t = trans[pos];
    t.src = pos;
    t.dest = pos + 1;
    t.lpTrans = 0;
    state[pos].out.push_back (pos);
    state[pos+1].in.push_back (pos);
    if (pos < seq.size()) {
      state[pos+1].lpAbsorb[seq[pos]] = 0;
      state[pos+1].alignPath[rowIndex].push_back (true);
    }
  }
}

Profile Profile::leftMultiply (gsl_matrix* sub) const {
  Profile prof (*this);
  for (ProfileStateIndex i = 0; i < size(); ++i)
    if (!state[i].isNull()) {
      for (AlphTok c = 0; c < alphSize; ++c) {
	LogProb lp = -numeric_limits<double>::infinity();
	for (AlphTok d = 0; d < alphSize; ++d)
	  lp = log_sum_exp (lp, gsl_matrix_get(sub,c,d) + state[i].lpAbsorb[d]);
	prof.state[i].lpAbsorb[c] = lp;
      }
    }
  return prof;
}

const ProfileTransition* Profile::getTrans (ProfileStateIndex src, ProfileStateIndex dest) const {
  for (auto t : state[dest].in)
    if (trans[t].src == src)
      return &trans[t];
  return NULL;
}
