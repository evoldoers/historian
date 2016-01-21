#include <iostream>
#include <iomanip>
#include "../src/alignpath.h"

int main (int argc, char **argv) {
  if (argc < 3) {
    cout << "Usage: " << argv[0] << " <align1> <align2> ...\n";
    exit (EXIT_FAILURE);
  }

  map<string,AlignRowIndex> nameToRowIndex;
  AlignRowIndex nextRowIndex = 0;
  vector<AlignPath> paths;

  for (int n = 2; n < argc; ++n) {
    vguard<FastSeq> fs = readFastSeqs (argv[n]);
    AlignPath tmpPath = gappedFastaToAlignPath (fs);
    AlignPath path;
    for (size_t n = 0; n < fs.size(); ++n) {
      if (nameToRowIndex.find(fs[n].name) == nameToRowIndex.end())
	nameToRowIndex[fs[n].name] = nextRowIndex++;
      path[nameToRowIndex[fs[n].name]] = tmpPath[n];
    }
    paths.push_back (path);
  }

  AlignPath merge = alignPathMerge (paths);
  for (AlignRowIndex r = 0; r < nextRowIndex; ++r) {
    cout << std::setw(4) << r << ' ';
    for (AlignColIndex c = 0; c < merge[r].size(); ++c)
      cout << (merge[r][c] ? 'x' : '-');
    cout << endl;
  }

  exit (EXIT_SUCCESS);
}
