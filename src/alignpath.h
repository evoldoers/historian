#ifndef ALIGNPATH_INCLUDED
#define ALIGNPATH_INCLUDED

#include <map>
#include <set>
#include "vguard.h"
#include "fastseq.h"

typedef size_t AlignRowIndex;
typedef size_t AlignColIndex;
typedef vguard<bool> AlignRowPath;
typedef map<AlignRowIndex,AlignRowPath> AlignPath;

AlignColIndex alignPathColumns (const AlignPath& a);
SeqIdx alignPathResiduesInRow (const AlignRowPath& r);

AlignPath alignPathUnion (const AlignPath& a1, const AlignPath& a2);  // simple union (no AlignRowIndex shared between a1 & a2)
AlignPath alignPathMerge (const vguard<AlignPath>& alignments);  // synchronized merge
AlignPath alignPathConcat (const AlignPath& a1, const AlignPath& a2);  // lengthwise concatenation
AlignPath alignPathConcat (const AlignPath& a1, const AlignPath& a2, const AlignPath& a3);

struct Alignment {
  vguard<FastSeq> ungapped;
  AlignPath path;
  Alignment (const vguard<FastSeq>& gapped);
  Alignment (const vguard<FastSeq>& ungapped, const AlignPath& path);
  vguard<FastSeq> gapped() const;
  static bool isGap (char c) { return c == '-' || c == '.'; }
};

struct GuideAlignmentDistanceMetric {
  // cumulativeMatches[row1][row2][col] = number of matches before column #col of pairwise alignment of (row1,row2)
  map<AlignRowIndex,map<AlignRowIndex,vguard<int> > > cumulativeMatches;
  // rowPosToCol[row][seqpos] = alignment column number of position #seqpos of row #row
  map<AlignRowIndex,vguard<AlignColIndex> > rowPosToCol;

  GuideAlignmentDistanceMetric (const AlignPath& guide);

  inline int pairwiseDistance (AlignRowIndex row1, SeqIdx pos1, AlignRowIndex row2, SeqIdx pos2) const {
    int minRow, maxRow;
    if (row1 < row2) { minRow = row1; maxRow = row2; }
    else { minRow = row2; maxRow = row1; }
    const vguard<int>& cm = cumulativeMatches.at(row1).at(row2);
    return abs (cm[rowPosToCol.at(row1).at(pos1)] - cm[rowPosToCol.at(row2).at(pos2)]);
  }

  inline int distance (const map<AlignRowIndex,SeqIdx>& coords1, const map<AlignRowIndex,SeqIdx>& coords2) const {
    int d = 0;
    for (auto& rp1 : coords1)
      for (auto& rp2 : coords2)
	d = max (d, pairwiseDistance (rp1.first, rp1.second, rp2.first, rp2.second));
    return d;
  }

  inline bool closeEnough (map<AlignRowIndex,SeqIdx>& coords1,
			   const map<AlignRowIndex,SeqIdx>& coords2,
			   int maxDistance) const {
    for (auto& rp1 : coords1)
      for (auto& rp2 : coords2)
	if (pairwiseDistance (rp1.first, rp1.second, rp2.first, rp2.second) > maxDistance)
	  return false;
    return true;
  }
};

#endif /* ALIGNPATH_INCLUDED */
