#ifndef ALIGNPATH_INCLUDED
#define ALIGNPATH_INCLUDED

#include <map>
#include <vector>

using namespace std;

typedef size_t AlignRowIndex;
typedef vector<bool> AlignRowPath;
typedef map<AlignRowIndex,AlignRowPath> AlignPath;

AlignPath unionAlignments (const AlignPath& a1, const AlignPath& a2);  // simple union (no AlignRowIndex shared between a1 & a2)
AlignPath syncAlignments (const vector<AlignPath>& alignments);  // synchronized merge
AlignPath concatAlignments (const AlignPath& a1, const AlignPath& a2);  // lengthwise concatenation

#endif /* ALIGNPATH_INCLUDED */
