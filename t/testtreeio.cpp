#include <iostream>
#include <fstream>
#include <string.h>
#include "../src/tree.h"
#include "../src/jsonutil.h"

int main (int argc, char **argv) {
  if (argc != 2 && argc != 3) {
    cout << "Usage: " << argv[0] << " <treefile> [<new root>]\n";
    exit (EXIT_FAILURE);
  }

  Tree tree;
  ifstream in (argv[1]);
  string s (JsonUtil::readStringFromStream (in));
  tree.parse (s);
  if (argc == 3)
    tree = tree.rerootAbove (string (argv[2]));
  cout << tree.toString() << endl;
  
  exit (EXIT_SUCCESS);
}
