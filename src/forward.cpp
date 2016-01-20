#include "forward.h"

ForwardMatrix::ForwardMatrix (const Profile& x, const Profile& y, const PairHMM& hmm)
  : x(x),
    y(y),
    hmm(hmm),
    alphSize (hmm.alphabetSize()),
    subx (x.leftMultiply (hmm.l.subMat)),
    suby (y.leftMultiply (hmm.r.subMat)),
    cellStorage (x.size() * y.size() * PairHMM::TotalStates, -numeric_limits<double>::infinity()),
    absorbScratch (hmm.alphabetSize()),
    insx (x.size(), -numeric_limits<double>::infinity()),
    insy (y.size(), -numeric_limits<double>::infinity())
{
  for (ProfileStateIndex i = 0; i < x.size(); ++i)
    if (!x.state[i].isNull())
      insx[i] = logInnerProduct (hmm.l.insVec, x.state[i].lpAbsorb);

  for (ProfileStateIndex j = 0; j < y.size(); ++j)
    if (!y.state[j].isNull())
      insy[j] = logInnerProduct (hmm.r.insVec, y.state[j].lpAbsorb);

  // WRITE ME: main DP loop
}

