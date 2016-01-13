#include <iostream>
#include <fstream>
#include <string.h>
#include "../src/logsumexp.h"

using namespace std;

int main (int argc, char **argv) {
  if (argc > 2 || (argc == 2 && strcmp(argv[1],"-slow") != 0 && strcmp(argv[1],"-fast") != 0)) {
    cout << "Usage: " << argv[0] << " [-slow|-fast]\n";
    exit (EXIT_FAILURE);
  }

  const bool slow = (argc == 2 && strcmp(argv[1],"-slow") == 0);
  cerr << "(running in " << (slow ? "slow" : "fast") << " mode)" << endl;
  for (double x = 0; x < 2; x += .1)
    for (double y = 0; y < 2; y += .1)
      cout << x << ' ' << y << ' ' << (slow ? log_sum_exp_slow(x,y) : log_sum_exp(x,y)) << endl;
  
  exit (EXIT_SUCCESS);
}
