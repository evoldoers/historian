#include "alignpath.h"

AlignPath unionAlignments (const AlignPath& a1, const AlignPath& a2) {
  AlignPath a = a1;
  a.insert (a2.begin(), a2.end());
  return a;
}

AlignPath concatAlignments (const AlignPath& a1, const AlignPath& a2) {
  AlignPath a = a1;
  for (auto& iter : a2) {
    const AlignRowIndex row = iter.first;
    const AlignRowPath& rPath = iter.second;
    AlignRowPath& lPath = a[row];
    lPath.insert (lPath.end(), rPath.begin(), rPath.end());
  }
  return a;
}
