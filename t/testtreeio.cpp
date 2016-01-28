#include <iostream>
#include <fstream>
#include <string.h>
#include "../src/tree.h"
#include "../src/jsonutil.h"

int main (int argc, char **argv) {
  if (argc != 2) {
    cout << "Usage: " << argv[0] << " <treefile>\n";
    exit (EXIT_FAILURE);
  }

  Tree tree;
  ifstream in (argv[1]);
  string s (JsonUtil::readStringFromStream (in));
  tree.parse (s);
  cout << tree.toString() << endl;
  
  exit (EXIT_SUCCESS);
}
