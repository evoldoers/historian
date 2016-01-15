#ifndef ALIGNPATH_INCLUDED
#define ALIGNPATH_INCLUDED

#include <map>
#include <set>
#include <vector>

typedef size_t AlignRowIndex;
typedef vector<bool> AlignRowPath;
typedef map<AlignRowIndex,AlignRowPath> AlignPath;

AlignPath combineAlignments (const set<AlignPath>& alignments);
AlignPath concatenateAlignments (const AlignPath& a1, const AlignPath& a2);

#endif /* ALIGNPATH_INCLUDED */
