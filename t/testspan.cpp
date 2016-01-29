#include <iostream>
#include <fstream>
#include <string.h>
#include "../src/jsonutil.h"
#include "../src/span.h"
#include "../src/logger.h"
#include "../src/diagenv.h"

int main (int argc, char **argv) {
  if (argc != 4) {
    cout << "Usage: " << argv[0] << " <seqfile> <modelfile> <time>\n";
    exit (EXIT_FAILURE);
  }

  const vguard<FastSeq> seqs = readFastSeqs (argv[1]);
  Require (seqs.size() >= 2, "Sequence file must have at least two sequences");

  RateModel rates;
  ifstream in (argv[2]);
  ParsedJson pj (in);
  rates.read (pj.value);

  const double time = atof (argv[3]);
  
  //  logger.setVerbose(6);

  ForwardMatrix::random_engine generator = ForwardMatrix::newRNG();
  DiagEnvParams dep;
  AlignGraph ag (seqs, rates, time, dep, generator);
  vguard<FastSeq> gapped = ag.mstGapped();
  writeFastaSeqs (cout, gapped);
  
  exit (EXIT_SUCCESS);
}
