#include <iostream>
#include <fstream>
#include <iomanip>
#include "../src/forward.h"
#include "../src/util.h"
#include "../src/jsonutil.h"

int main (int argc, char **argv) {
  if (argc != 4 && argc != 5) {
    cout << "Usage: " << argv[0] << " <sequences> <modelfile> <xtime> [<ytime>]\n";
    exit (EXIT_FAILURE);
  }

  vguard<FastSeq> seqs = readFastSeqs (argv[1]);
  Assert (seqs.size() == 2, "Expected two sequences in file %s", argv[1]);

  RateModel rates;
  ifstream in (argv[2]);
  ParsedJson pj (in);
  rates.read (pj.value);

  ProbModel xprobs (rates, atof (argv[3]));
  ProbModel yprobs (rates, atof (argv[argc > 4 ? 4 : 3]));
  gsl_vector* eqm = rates.getEqmProb();
  PairHMM hmm (xprobs, yprobs, eqm);

  Profile xprof (rates.alphabet, seqs[0], 1);
  Profile yprof (rates.alphabet, seqs[1], 2);
  ForwardMatrix forward (xprof, yprof, hmm, 0);

  set<ForwardMatrix::CellCoords> allCells;
  allCells.insert (forward.startCell);
  allCells.insert (forward.endCell);
  for (ProfileStateIndex xpos = 0; xpos < xprof.size() - 1; ++xpos)
    for (ProfileStateIndex ypos = 0; ypos < yprof.size() - 1; ++ypos)
      for (PairHMM::State s : hmm.states())
	if (xpos > 0 || ypos > 0)
	  allCells.insert (ForwardMatrix::CellCoords (xpos, ypos, s));

  Profile prof = forward.makeProfile (allCells, true);
  prof.calcSumPathAbsorbProbs (hmm.logRoot);
  prof.writeJson (cout);

  exit (EXIT_SUCCESS);
}
