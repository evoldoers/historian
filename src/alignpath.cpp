#include "alignpath.h"

AlignPath mergeAlignments (const AlignPath& a1, const AlignPath& a2) {
  AlignPath a = a1;
  a.insert (a2.begin(), a2.end());
  return a;
}

