#include <iostream>
#include <fstream>
#include <string.h>
#include "../src/profile.h"

int main (int argc, char **argv) {
  if (argc != 3) {
    cout << "Usage: " << argv[0] << " <alphabet> <sequence>\n";
    exit (EXIT_FAILURE);
  }

  FastSeq fs;
  
  const string alphabet (argv[1]);
  fs.seq = (argv[2]);
  const vguard<AlphTok> dsq = fs.tokens (alphabet);

  Profile prof (alphabet.size(), dsq, 0);
  prof.writeJson (cout);
  
  exit (EXIT_SUCCESS);
}
