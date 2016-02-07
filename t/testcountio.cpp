#include <iostream>
#include <fstream>
#include <string.h>
#include "../src/model.h"
#include "../src/jsonutil.h"

int main (int argc, char **argv) {
  if (argc != 2) {
    cout << "Usage: " << argv[0] << " <countfile>\n";
    exit (EXIT_FAILURE);
  }

  ifstream in (argv[1]);
  ParsedJson pj (in);
  EventCounts c;
  c.read (pj.value);
  c.writeJson (cout);
  
  exit (EXIT_SUCCESS);
}
