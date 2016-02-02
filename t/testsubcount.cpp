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
  gsl_matrix *count = gsl_matrix_calloc (rates.alphabetSize(), rates.alphabetSize());

  eigen.accumSubCounts (count, src, dest, 1, sub, esub);

  gsl_matrix_free (sub);
  gsl_matrix_complex_free (esub);

  gsl_vector *root = gsl_vector_calloc (rates.alphabetSize());
  gsl_vector_set (root, src, 1);

  rates.writeSubCounts (cout, root, count);

  gsl_matrix_free (count);
  gsl_vector_free (root);

  exit (EXIT_SUCCESS);
}
