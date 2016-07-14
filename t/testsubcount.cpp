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

  vguard<gsl_matrix*> sub = eigen.getSubProbMatrix(t);
  vguard<gsl_matrix_complex*> esub = eigen.eigenSubCount(t);
  vguard<vguard<vguard<double> > > count (1, vguard<vguard<double> > (rates.alphabetSize(), vguard<double> (rates.alphabetSize(), 0)));

  eigen.accumSubCounts (0, count[0], src, dest, 1, sub[0], esub[0]);

  for (auto& s: sub)
    gsl_matrix_free (s);
  for (auto& e: esub)
    gsl_matrix_complex_free (e);

  vguard<vguard<double> > root (1, vguard<double> (rates.alphabetSize(), 0));
  root[0][src] = 1;

  rates.writeSubCounts (cout, root, count);
  cout << endl;

  exit (EXIT_SUCCESS);
}
