#include <iostream>
#include <fstream>
#include <string.h>
#include "../src/nexus.h"

int main (int argc, char **argv) {
  if (argc != 2) {
    cout << "Usage: " << argv[0] << " <nexus file>\n";
    exit (EXIT_FAILURE);
  }

  ifstream in (argv[1]);
  const NexusData nexus (in);
  nexus.write (cout);
  
  exit (EXIT_SUCCESS);
}
