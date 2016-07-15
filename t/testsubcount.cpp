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
  vguard<vguard<vguard<double> > > count (rates.components(), vguard<vguard<double> > (rates.alphabetSize(), vguard<double> (rates.alphabetSize(), 0)));

  vguard<double> p (rates.components());
  double norm = 0;
  for (int cpt = 0; cpt < rates.components(); ++cpt)
    norm += (p[cpt] = gsl_matrix_get (sub[cpt], src, dest));
  for (int cpt = 0; cpt < rates.components(); ++cpt)
    eigen.accumSubCounts (cpt, count[cpt], src, dest, p[cpt] / norm, sub[cpt], esub[cpt]);

  for (auto& s: sub)
    gsl_matrix_free (s);
  for (auto& e: esub)
    gsl_matrix_complex_free (e);

  vguard<vguard<double> > root (rates.components(), vguard<double> (rates.alphabetSize(), 0));
  for (int cpt = 0; cpt < rates.components(); ++cpt)
    root[cpt][src] = p[cpt] / norm;

  rates.writeSubCounts (cout, root, count);
  cout << endl;

  exit (EXIT_SUCCESS);
}
