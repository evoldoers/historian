#include <iostream>
#include <fstream>
#include <string.h>
#include "../src/model.h"
#include "../src/jsonutil.h"
#include "../src/sumprod.h"
#include "../src/logger.h"

int main (int argc, char **argv) {
  if (argc != 5) {
    cout << "Usage: " << argv[0] << " <modelfile> <src> <dest> <time>\n";
    exit (EXIT_FAILURE);
  }

  RateModel rates;
  ifstream in (argv[1]);
  ParsedJson pj (in);
  rates.read (pj.value);

  const AlphTok src = rates.tokenize (argv[2][0]);
  const AlphTok dest = rates.tokenize (argv[3][0]);
  const double t = atof (argv[4]);

  //  logger.setVerbose (8);
  
  EigenModel eigen (rates);

  gsl_matrix *sub = eigen.getSubProbMatrix(t);
  gsl_matrix_complex *esub = eigen.eigenSubCount(t);
  vguard<vguard<double> > count (rates.alphabetSize(), vguard<double> (rates.alphabetSize(), 0));

  eigen.accumSubCounts (count, src, dest, 1, sub, esub);

  gsl_matrix_free (sub);
  gsl_matrix_complex_free (esub);

  vguard<double> root (rates.alphabetSize(), 0);
  root[src] = 1;

  rates.writeSubCounts (cout, root, count);

  exit (EXIT_SUCCESS);
}
