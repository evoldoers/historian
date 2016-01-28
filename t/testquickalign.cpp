#include <iostream>
#include <fstream>
#include <string.h>
#include "../src/jsonutil.h"
#include "../src/quickalign.h"
#include "../src/logger.h"

int main (int argc, char **argv) {
  if (argc != 4) {
    cout << "Usage: " << argv[0] << " <seqfile> <modelfile> <time>\n";
    exit (EXIT_FAILURE);
  }

  const vguard<FastSeq> seqs = readFastSeqs (argv[1]);
  Require (seqs.size() == 2, "Sequence file must have exactly two sequences");

  RateModel rates;
  ifstream in (argv[2]);
  ParsedJson pj (in);
  rates.read (pj.value);

  const double time = atof (argv[3]);
  
  //  logger.setVerbose(6);

  DiagonalEnvelope env (seqs[0], seqs[1]);
  env.initFull();

  QuickAlignMatrix mx (env, rates, time);
  Alignment align = mx.alignment();
  vguard<FastSeq> gapped = align.gapped();
  writeFastaSeqs (cout, gapped);
  
  exit (EXIT_SUCCESS);
}
