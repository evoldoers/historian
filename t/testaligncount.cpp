#include <iostream>
#include <fstream>
#include <string.h>
#include "../src/model.h"
#include "../src/jsonutil.h"
#include "../src/sumprod.h"
#include "../src/logger.h"

int main (int argc, char **argv) {
  if (argc != 4) {
    cout << "Usage: " << argv[0] << " <model> <alignment> <tree>\n";
    exit (EXIT_FAILURE);
  }

  RateModel rates;
  ifstream in (argv[1]);
  ParsedJson pj (in);
  rates.read (pj.value);

  const vguard<FastSeq> gapped = readFastSeqs (argv[2]);

  ifstream treeStream (argv[3]);
  Tree tree (JsonUtil::readStringFromStream (treeStream));

  //  logger.setVerbose (8);
  
  AlignColSumProduct colSumProd (rates, tree, gapped);

  gsl_matrix_complex *eigenCount = gsl_matrix_complex_calloc (rates.alphabetSize(), rates.alphabetSize());
  gsl_vector *root = gsl_vector_calloc (rates.alphabetSize());

  while (!colSumProd.alignmentDone()) {
    colSumProd.fillUp();
    colSumProd.fillDown();
    colSumProd.accumulateCounts (root, eigenCount);
    colSumProd.nextColumn();
  }

  gsl_matrix *count = colSumProd.getSubCounts (eigenCount);

  rates.writeSubCounts (cout, root, count);

  gsl_matrix_complex_free (eigenCount);
  gsl_matrix_free (count);
  gsl_vector_free (root);

  exit (EXIT_SUCCESS);
}
