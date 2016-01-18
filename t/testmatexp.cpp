#include <iostream>
#include <fstream>
#include <string.h>
#include "../src/model.h"
#include "../src/jsonutil.h"

int main (int argc, char **argv) {
  if (argc != 3) {
    cout << "Usage: " << argv[0] << " <modelfile> <time>\n";
    exit (EXIT_FAILURE);
  }

  RateModel rates;
  ifstream in (argv[1]);
  ParsedJson pj (in);
  rates.read (pj.value);

  ProbModel probs (rates, atof (argv[2]));
  probs.write (cout);
  
  exit (EXIT_SUCCESS);
}
