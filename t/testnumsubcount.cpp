#include <iostream>
#include <fstream>
#include <string.h>
#include "../src/model.h"
#include "../src/jsonutil.h"
#include "../src/sumprod.h"
#include "../src/logger.h"

int main (int argc, char **argv) {
  if (argc != 7) {
    cout << "Usage: " << argv[0] << " <modelfile> <a> <b> <i> <j> <time>\n";
    exit (EXIT_FAILURE);
  }

  RateModel rates;
  ifstream in (argv[1]);
  ParsedJson pj (in);
  rates.read (pj.value);

  const AlphTok a = rates.tokenize (argv[2][0]);
  const AlphTok b = rates.tokenize (argv[3][0]);
  const AlphTok i = rates.tokenize (argv[4][0]);
  const AlphTok j = rates.tokenize (argv[5][0]);
  const double T = atof (argv[6]);

  logger.setVerbose (8);
  
  EigenModel eigen (rates);

  gsl_matrix *sub = eigen.getSubProbMatrix(T);
  gsl_matrix_complex *esub = eigen.eigenSubCount(T);

  const double count = eigen.getSubCount (a, b, i, j, sub, esub);

  const double nSteps = 1e6;
  const double tStep = 1 / nSteps;
  double numCount = 0;
  for (double t = tStep; t < T; t += tStep)
    numCount += eigen.getSubProb(t,a,i) * eigen.getSubProb(T-t-tStep,j,b);
  numCount *= gsl_matrix_get (sub, i, j) * tStep / eigen.getSubProb(T,a,b);

  cout << "Eigenvector method: " << count << endl;
  cout << "Numerical integration: " << numCount << endl;
  
  gsl_matrix_free (sub);
  gsl_matrix_complex_free (esub);
  
  exit (EXIT_SUCCESS);
}
