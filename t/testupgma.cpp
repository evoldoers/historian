#include <iostream>
#include <fstream>
#include <string.h>
#include "../src/model.h"
#include "../src/jsonutil.h"
#include "../src/tree.h"
#include "../src/logger.h"

int main (int argc, char **argv) {
  if (argc != 3) {
    cout << "Usage: " << argv[0] << " <modelfile> <alignfile>\n";
    exit (EXIT_FAILURE);
  }

  RateModel rates;
  ifstream in (argv[1]);
  ParsedJson pj (in);
  rates.read (pj.value);

  const vguard<FastSeq> gapped = readFastSeqs (argv[2]);

  //  logger.setVerbose(6);
  auto dist = rates.distanceMatrix (gapped);
  Tree tree;
  tree.buildByUPGMA (gapped, dist);

  cout << tree.toString() << endl;
  
  exit (EXIT_SUCCESS);
}
