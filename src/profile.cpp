#include "profile.h"

ProfileTransition::ProfileTransition (AlphTok alphSize)
  : src (0),
    dest (0),
    absorb (alphSize, -numeric_limits<double>::infinity())
{ }

Profile::Profile (AlphTok alphSize, const vguard<AlphTok>& seq)
  : alphSize (alphSize),
    state (seq.size() + 1),
    trans (seq.size(), ProfileTransition (alphSize))
{
  for (size_t pos = 0; pos < seq.size(); ++pos) {
    trans[pos].src = pos;
    trans[pos].dest = pos + 1;
    trans[pos].absorb[seq[pos]] = 0;
    state[pos].out.push_back (pos);
    state[pos+1].in.push_back (pos);
  }
}
