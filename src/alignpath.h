#ifndef ALIGNPATH_INCLUDED
#define ALIGNPATH_INCLUDED

#include <map>
#include <set>
#include <vector>
#include "fastseq.h"

typedef size_t AlignRowIndex;
typedef size_t AlignColIndex;
typedef vector<bool> AlignRowPath;
typedef map<AlignRowIndex,AlignRowPath> AlignPath;

AlignPath alignPathUnion (const AlignPath& a1, const AlignPath& a2);  // simple union (no AlignRowIndex shared between a1 & a2)
AlignPath alignPathMerge (const vector<AlignPath>& alignments);  // synchronized merge
AlignPath alignPathConcat (const AlignPath& a1, const AlignPath& a2);  // lengthwise concatenation
AlignPath alignPathConcat (const AlignPath& a1, const AlignPath& a2, const AlignPath& a3);

struct Alignment {
  vector<FastSeq> ungapped;
  AlignPath path;
  Alignment (const vector<FastSeq>& gapped);
  Alignment (const vector<FastSeq>& ungapped, const AlignPath& path);
  vector<FastSeq> gapped() const;
  static bool isGap (char c) { return c == '-' || c == '.'; }
};

#endif /* ALIGNPATH_INCLUDED */
