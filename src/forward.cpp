#include "forward.h"

ForwardMatrix::ForwardMatrix (const Profile& x, const Profile& y, const PairHMM& hmm)
  : x(x),
    y(y),
    hmm(hmm),
    alphSize (hmm.alphabetSize()),
    xSize (x.size()),
    ySize (y.size()),
    subx (x.leftMultiply (hmm.l.subMat)),
    suby (y.leftMultiply (hmm.r.subMat)),
    cellStorage (x.size() * y.size() * PairHMM::TotalStates, -numeric_limits<double>::infinity()),
    absorbScratch (hmm.alphabetSize()),
    insx (x.size(), -numeric_limits<double>::infinity()),
    insy (y.size(), -numeric_limits<double>::infinity())
{
  for (ProfileStateIndex i = 1; i < xSize - 1; ++i)
    insx[i] = logInnerProduct (hmm.l.insVec, x.state[i].lpAbsorb);

  for (ProfileStateIndex j = 1; j < ySize - 1; ++j)
    insy[j] = logInnerProduct (hmm.r.insVec, y.state[j].lpAbsorb);

  cell(0,0,IMM) = 0;
  for (ProfileStateIndex i = 0; i < xSize - 1; ++i)
    for (ProfileStateIndex j = 0; j < ySize - 1; ++j) {
      if (i > 0) {
	// WRITE ME: transitions into IMD, IIW, IIX
      }

      if (j > 0) {
	// WRITE ME: transitions into IDM, IMI, IDI
      }

      if (i > 0 && j > 0) {
	// WRITE ME: transitions into IMM
      }
    }

  // WRITE ME: transitions into EEE
}

