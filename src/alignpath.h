#ifndef ALIGNPATH_INCLUDED
#define ALIGNPATH_INCLUDED

#include <map>
#include <set>
#include <vector>

using namespace std;

typedef size_t AlignRowIndex;
typedef vector<bool> AlignRowPath;
typedef map<AlignRowIndex,AlignRowPath> AlignPath;

AlignPath mergeAlignments (const AlignPath& a1, const AlignPath& a2);  // simple union (no AlignRowIndex shared between a1 & a2)
AlignPath combineAlignments (const set<AlignPath>& alignments);  // synchronized mesh
AlignPath concatenateAlignments (const AlignPath& a1, const AlignPath& a2);  // lengthwise concatenation

#endif /* ALIGNPATH_INCLUDED */
