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


// sketch of syncAlignments algorithm:
// build bidirectional map from (align#,column#) <==> (row#,residue#)
// define function linkedColumns: maps single (align#,column#) to set of { (align#,column#) }
// create cursor: a vector of nextColumn#'s indexed by align#
// define readyToAdvance(align#,cursor): linkedColumns are all cursor-ready
