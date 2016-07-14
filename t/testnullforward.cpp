#include <iostream>
#include <fstream>
#include <iomanip>
#include "../src/forward.h"
#include "../src/util.h"
#include "../src/jsonutil.h"

int main (int argc, char **argv) {

  if (argc != 3 && argc != 4) {
    cout << "Usage: " << argv[0] << " <modelfile> <xtime> [<ytime>]\n";
    exit (EXIT_FAILURE);
  }

  RateModel rates;
  ifstream in (argv[1]);
  ParsedJson pj (in);
  rates.read (pj.value);

  ProbModel xprobs (rates, atof (argv[2]));
  ProbModel yprobs (rates, atof (argv[argc > 3 ? 3 : 2]));
  vguard<gsl_vector*> eqm = rates.insProb;
  PairHMM hmm (xprobs, yprobs, eqm);

  FastSeq x, y;
  x.name = "x";
  x.seq = "acg";

  y.name = "y";
  y.seq = "cag";
  
  Profile xprof (1, rates.alphabet, x, 1);
  Profile yprof (1, rates.alphabet, y, 2);

  xprof.state[2].lpAbsorb.clear();

  yprof.state[1].lpAbsorb.clear();

  ForwardMatrix forward (xprof, yprof, hmm, 0, GuideAlignmentEnvelope());

  set<ForwardMatrix::CellCoords> allCells;
  allCells.insert (forward.startCell);
  allCells.insert (forward.endCell);
  for (ProfileStateIndex xpos = 0; xpos < xprof.size() - 1; ++xpos)
    for (ProfileStateIndex ypos = 0; ypos < yprof.size() - 1; ++ypos)
      for (PairHMM::State s : hmm.states())
	if (xpos > 0 || ypos > 0)
	  allCells.insert (ForwardMatrix::CellCoords (xpos, ypos, s));

  Profile prof = forward.makeProfile (allCells, ForwardMatrix::KeepAll);
  prof.calcSumPathAbsorbProbs (vguard<LogProb>(1,0), hmm.logRoot);
  prof.writeJson (cout);

  exit (EXIT_SUCCESS);
}
