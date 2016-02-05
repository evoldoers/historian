#include <iostream>
#include <fstream>
#include <string.h>
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

  tree.reorder (gapped);
  
  //  logger.setVerbose (8);
  
  AlignColSumProduct colSumProd (rates, tree, gapped);

  gsl_matrix_complex *eigenCount = gsl_matrix_complex_calloc (rates.alphabetSize(), rates.alphabetSize());
  vguard<vguard<double> > count (rates.alphabetSize(), vguard<double> (rates.alphabetSize(), 0));
  vguard<double> root (rates.alphabetSize(), 0);

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
    count = colSumProd.getSubCounts (eigenCount);
  rates.writeSubCounts (cout, root, count);

  gsl_matrix_complex_free (eigenCount);

  exit (EXIT_SUCCESS);
}
