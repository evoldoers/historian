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
  vguard<gsl_vector*> eqm = rates.insProb;
  PairHMM hmm (xprobs, yprobs, eqm);

  Profile xprof (1, rates.alphabet, seqs[0], 1);
  Profile yprof (1, rates.alphabet, seqs[1], 2);
  ForwardMatrix forward (xprof, yprof, hmm, 0, GuideAlignmentEnvelope());
  BackwardMatrix backward (forward);

  cout << "Forward score: " << forward.lpEnd << endl;
  cout << "Backward score: " << backward.lpStart() << endl;

  auto bestCells = backward.cellsAbovePostProbThreshold (.5);
  while (!bestCells.empty()) {
    cout << "P" << backward.cellName (bestCells.top()) << " = " << exp(bestCells.top().logPostProb) << endl;
    bestCells.pop();
  }
  
  exit (EXIT_SUCCESS);
}
