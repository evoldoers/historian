#include <iostream>
#include <fstream>
#include <string.h>
#include "../src/fastseq.h"

int main (int argc, char **argv) {
  if (argc != 2) {
    cout << "Usage: " << argv[0] << " <seqfile>\n";
    exit (EXIT_FAILURE);
  }

  const vguard<FastSeq> gapped = readFastSeqs (argv[1]);
  writeFastaSeqs (cout, gapped);

  exit (EXIT_SUCCESS);
}
