#include <iostream>
#include <iomanip>
#include "../src/alignpath.h"

int main (int argc, char **argv) {
  if (argc < 3) {
    cout << "Usage: " << argv[0] << " <align1> <align2> ...\n";
    exit (EXIT_FAILURE);
  }

  map<string,AlignRowIndex> nameToRowIndex;
  vguard<FastSeq> ungapped;
  vguard<AlignPath> paths;

  for (int n = 1; n < argc; ++n) {
    vguard<FastSeq> gapped = readFastSeqs (argv[n]);
    Alignment align (gapped);
    AlignPath path;
    for (size_t n = 0; n < gapped.size(); ++n) {
      if (nameToRowIndex.find(gapped[n].name) == nameToRowIndex.end()) {
	nameToRowIndex[gapped[n].name] = ungapped.size();
	ungapped.push_back (align.ungapped[n]);
      }
      path[nameToRowIndex[gapped[n].name]] = align.path[n];
    }
    paths.push_back (path);
  }

  const AlignPath path = alignPathMerge (paths);
  const Alignment align (ungapped, path);
  writeFastaSeqs (cout, align.gapped());

  exit (EXIT_SUCCESS);
}
