#include <iostream>
#include <fstream>
#include <string.h>
#include <gsl/gsl_complex_math.h>
#include "../src/model.h"
#include "../src/jsonutil.h"
#include "../src/sumprod.h"
#include "../src/logger.h"

int main (int argc, char **argv) {
  if (argc != 5) {
    cout << "Usage: " << argv[0] << " [-eigen|-sub] <model> <alignment> <tree>\n";
    exit (EXIT_FAILURE);
  }

  const bool useEigen = strcmp(argv[1],"-eigen") == 0;
  Assert (useEigen || strcmp(argv[1],"-sub") == 0, "First argument must be -sub or -eigen");
  
  RateModel rates;
  ifstream in (argv[2]);
  ParsedJson pj (in);
  rates.read (pj.value);

  vguard<FastSeq> gapped = readFastSeqs (argv[3]);

  ifstream treeStream (argv[4]);
  Tree tree (JsonUtil::readStringFromStream (treeStream));

  tree.reorderSeqs (gapped);
  
  //  logger.setVerbose (8);
  
  AlignColSumProduct colSumProd (rates, tree, gapped);

  vguard<vguard<vguard<gsl_complex> > > eigenCount (1, vguard<vguard<gsl_complex> > (rates.alphabetSize(), vguard<gsl_complex> (rates.alphabetSize(), gsl_complex_rect(0,0))));
  vguard<vguard<vguard<double> > > count (1, vguard<vguard<double> > (rates.alphabetSize(), vguard<double> (rates.alphabetSize(), 0)));
  vguard<vguard<double> > root (1, vguard<double> (rates.alphabetSize(), 0));

  while (!colSumProd.alignmentDone()) {
    colSumProd.fillUp();
    colSumProd.fillDown();
    if (useEigen)
      colSumProd.accumulateEigenCounts (root, eigenCount);
    else
      colSumProd.accumulateSubCounts (root, count);
    colSumProd.nextColumn();
  }

  if (useEigen)
    count = colSumProd.eigen.getSubCounts (eigenCount);

  rates.writeSubCounts (cout, root, count);
  cout << endl;
  
  exit (EXIT_SUCCESS);
}
