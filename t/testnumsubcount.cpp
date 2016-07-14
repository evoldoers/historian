#include <iostream>
#include <fstream>
#include <string.h>
#include "../src/model.h"
#include "../src/jsonutil.h"
#include "../src/sumprod.h"
#include "../src/logger.h"

double jcProb (double lambda, double t, int i, int j) {
  const double e = exp(-lambda*t);
  return i == j ? (e + (1 - e) / 4) : ((1 - e) / 4);
}

int main (int argc, char **argv) {
  if (argc != 7 && argc != 8) {
    cout << "Usage: " << argv[0] << " <modelfile> <a> <b> <i> <j> <time> [<lambda>]\n";
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

  //  logger.setVerbose (8);
  
  EigenModel eigen (rates);

  vguard<gsl_matrix*> sub = eigen.getSubProbMatrix(T);
  vguard<gsl_matrix_complex*> esub = eigen.eigenSubCount(T);

  const double count = eigen.getSubCount (0, a, b, i, j, sub[0], esub[0]);

  const double nSteps = 1e5;
  const double tStep = T / nSteps;
  double numCount = 0;
  for (double t = 0; t < T; t += tStep)
    numCount += eigen.getSubProb(0,t,a,i) * eigen.getSubProb(0,T-t-tStep,j,b);
  numCount *= gsl_matrix_get (rates.subRate[0], i, j) * tStep / eigen.getSubProb(0,T,a,b);

  cout << "Eigenvector method: " << count << endl;
  cout << "Numerical integration: " << numCount << endl;

  if (argc == 8) {
    Assert (i != j, "Need i!=j for exact Jukes-Cantor test");
    
    const double lambda = atof (argv[7]);
    if (a != i && j != b && a != b) {
      const double jcCount = (lambda/16) * (T + (2/lambda)*(exp(-lambda*T)-1) + T*exp(-lambda*T)) / (1 - exp(-lambda*T));
  
      cout << "Jukes-Cantor (lambda=" << lambda << "): " << jcCount << endl;
    }

    double jcNumCount = 0;
    for (double t = 0; t < T; t += tStep)
      jcNumCount += jcProb(lambda,t,a,i) * (lambda/4) * tStep * jcProb(lambda,T-t,j,b);
    jcNumCount /= jcProb(lambda,T,a,b);
    cout << "Jukes-Cantor numerical (lambda=" << lambda << "): " << jcNumCount << endl;

    cout << "Rate(i->j): " << gsl_matrix_get (rates.subRate[0], i, j) << endl;

    cout << "Eigen: P(a->i|T/3): " << eigen.getSubProb(0,T/3,a,i) << endl;
    cout << "Eigen: P(j->b|2T/3): " << eigen.getSubProb(0,2*T/3,j,b) << endl;
    cout << "Eigen: P(a->b|T): " << eigen.getSubProb(0,T,a,b) << endl;
    
    cout << "JC exact: P(a->i|T/3): " << jcProb(lambda,T/3,a,i) << endl;
    cout << "JC exact: P(j->b|2T/3): " << jcProb(lambda,2*T/3,j,b) << endl;
    cout << "JC exact: P(a->b|T): " << jcProb(lambda,T,a,b) << endl;
  }

  for (auto& s: sub)
    gsl_matrix_free (s);
  for (auto& e: esub)
    gsl_matrix_complex_free (e);
  
  exit (EXIT_SUCCESS);
}
