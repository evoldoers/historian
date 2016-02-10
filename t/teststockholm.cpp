#include <iostream>
#include <fstream>
#include <string.h>
#include "../src/stockholm.h"

int main (int argc, char **argv) {
  if (argc != 2) {
    cout << "Usage: " << argv[0] << " <Stockholm file>\n";
    exit (EXIT_FAILURE);
  }

  ifstream in (argv[1]);
  const Stockholm stock (in);
  stock.write (cout);
  
  exit (EXIT_SUCCESS);
}
