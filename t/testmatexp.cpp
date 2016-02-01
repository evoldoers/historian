#include <iostream>
#include <fstream>
#include <string.h>
#include "../src/model.h"
#include "../src/jsonutil.h"
#include "../src/sumprod.h"
#include "../src/logger.h"

int main (int argc, char **argv) {
  const bool useEigen = argc == 4 && string(argv[1]) == "-eigen";
  if (argc != 3 && !useEigen) {
    cout << "Usage: " << argv[0] << " [-eigen] <modelfile> <time>\n";
    exit (EXIT_FAILURE);
  }

  RateModel rates;
  ifstream in (argv[argc-2]);
  ParsedJson pj (in);
  rates.read (pj.value);

  const double t = atof (argv[argc-1]);

  ProbModel probs (rates, t);
  if (useEigen) {
    //    logger.setVerbose(8);
    cerr << "Using eigendecomposition to compute matrix exponential" << endl;
    EigenModel eigen (rates);
    gsl_matrix_free (probs.subMat);
    probs.subMat = eigen.getSubProb (t);
  }
  probs.write (cout);
  
  exit (EXIT_SUCCESS);
}
