#include "alignpath.h"
#include "util.h"

// map used by alignPathMerge
struct AlignSeqMap {
  typedef size_t AlignNum;
  const vector<AlignPath>& alignments;
  map<AlignRowIndex,SeqIdx> seqLen;
  vector<AlignColIndex> alignCols;
  map<AlignNum,map<AlignColIndex,map<AlignRowIndex,SeqIdx> > > alignColRowToPos;
  map<AlignRowIndex,map<SeqIdx,map<AlignNum,AlignColIndex> > > rowPosAlignToCol;
  AlignSeqMap (const vector<AlignPath>& alignments);
  map<AlignNum,AlignColIndex> linkedColumns (AlignNum nAlign, AlignColIndex col) const;
};

AlignPath alignPathUnion (const AlignPath& a1, const AlignPath& a2) {
  AlignPath a = a1;
  a.insert (a2.begin(), a2.end());
  return a;
}

AlignPath alignPathConcat (const AlignPath& a1, const AlignPath& a2) {
  AlignPath a = a1;
  for (auto& iter : a1)
    Assert (a2.find(iter.first) != a2.end(), "Alignment row set mismatch");
  for (auto& iter : a2) {
    Assert (a1.find(iter.first) != a1.end(), "Alignment row set mismatch");
    const AlignRowIndex row = iter.first;
    const AlignRowPath& rPath = iter.second;
    AlignRowPath& lPath = a[row];
    lPath.insert (lPath.end(), rPath.begin(), rPath.end());
  }
  return a;
}

AlignSeqMap::AlignSeqMap (const vector<AlignPath>& alignments)
  : alignments (alignments)
{
  // get row indices and sequence lengths; confirm row & sequence lengths match
  for (auto& align : alignments) {
    if (align.size() == 0)
      alignCols.push_back (0);
    else {
      AlignColIndex* cols = NULL;
      for (auto& row_path : align) {
	const AlignRowIndex row = row_path.first;
	const AlignRowPath& path = row_path.second;
	if (cols == NULL) {
	  alignCols.push_back (path.size());
	  cols = &alignCols.back();
	} else
	  Assert (*cols == path.size(), "Incompatible number of columns in row #%d of alignment (%d != %d)", row, *cols, path.size());
	SeqIdx len = 0;
	for (bool b : path)
	  if (b)
	    ++len;
	if (seqLen.find(row) == seqLen.end())
	  seqLen[row] = len;
	else
	  Assert (seqLen[row] == len, "Incompatible number of residues for row #%d of alignment (%d != %d)", row, seqLen[row], len);
      }
    }
  }

  // build bidirectional map from (align#,column#) <==> (row#,residue#)
  for (size_t nAlign = 0; nAlign < alignments.size(); ++nAlign) {
    auto& align = alignments[nAlign];
    map<AlignRowIndex,SeqIdx> rowPos;
    for (auto& row_path : align)
      rowPos[row_path.first] = 0;
    for (AlignColIndex col = 0; col < alignCols[nAlign]; ++col) {
      rowPos.clear();
      for (auto& row_path : align)
	if (row_path.second[col]) {
	  const SeqIdx pos = rowPos[row_path.first]++;
	  alignColRowToPos[nAlign][col][row_path.first] = pos;
	  rowPosAlignToCol[row_path.first][pos][nAlign] = col;
	}
    }
  }
}

map<AlignSeqMap::AlignNum,AlignColIndex> AlignSeqMap::linkedColumns (AlignNum nAlign, AlignColIndex col) const {
  map<AlignNum,AlignColIndex> ac;
  for (auto& row_pos : alignColRowToPos.at(nAlign).at(col))
    for (auto& nAlign_col : rowPosAlignToCol.at(row_pos.first).at(row_pos.second)) {
      if (ac.find (nAlign_col.first) != ac.end())
	Assert (ac[nAlign_col.first] == nAlign_col.second, "Inconsistent alignments");
      ac.insert (nAlign_col);
    }
  return ac;
}

AlignPath alignPathMerge (const vector<AlignPath>& alignments) {
  const AlignSeqMap alignSeqMap (alignments);
  AlignPath a;
  for (auto& row_seqlen : alignSeqMap.seqLen)
    a[row_seqlen.first].clear();
  vector<AlignColIndex> nextCol (alignments.size(), 0);
  bool allDone, noneReady;
  do {
    allDone = noneReady = true;
    map<AlignSeqMap::AlignNum,AlignColIndex> linkedCols;
    for (AlignSeqMap::AlignNum n = 0; n < alignments.size(); ++n)
      if (nextCol[n] < alignSeqMap.alignCols[n]) {
	allDone = false;
	bool ready = true;
	linkedCols = alignSeqMap.linkedColumns (n, nextCol[n]);
	for (const auto& nAlign_col : linkedCols)
	  if (nextCol[nAlign_col.first] != nAlign_col.second) {
	    ready = false;
	    break;
	  }
	if (ready) {
	  noneReady = false;
	  if (linkedCols.size()) {
	    for (auto& idx_path : a)
	      idx_path.second.push_back (false);
	    for (const auto& nAlign_col : linkedCols)
	      for (const auto& row_path : alignments.at(nAlign_col.first))
		if (alignments.at(nAlign_col.first).at(row_path.first).at(nAlign_col.second))
		  a[row_path.first].back() = true;
	  } else
	    ++nextCol[n];  // empty column
	  break;
	}
      }
    Assert (!noneReady, "%s fail", __func__);
  } while (!allDone);
  return a;
}

AlignPath gappedFastaToAlignPath (const vector<FastSeq>& fs) {
  AlignPath align;
  for (AlignRowIndex row = 0; row < fs.size(); ++row) {
    AlignRowPath rowPath (fs[row].length());
    for (AlignColIndex col = 0; col < rowPath.size(); ++col)
      rowPath[col] = fs[row].seq[col] != '-' && fs[row].seq[col] != '.';
    align[row] = rowPath;
  }
  return align;
}
