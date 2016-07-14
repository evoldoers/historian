#include <iostream>
#include <fstream>
#include <string.h>
#include "../src/profile.h"

int main (int argc, char **argv) {
  if (argc != 3) {
    cout << "Usage: " << argv[0] << " <alphabet> <sequence>\n";
    exit (EXIT_FAILURE);
  }

  const string alphabet (argv[1]);
  FastSeq fs;
  fs.seq = (argv[2]);

  Profile prof (1, alphabet, fs, 0);
  prof.writeJson (cout);
  
  exit (EXIT_SUCCESS);
}
