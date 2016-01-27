#include <set>
#include <algorithm>
#include <iostream>
#include <cmath>
#include "diagenv.h"
#include "logger.h"

// require at least this ratio of sequenceLength/(kmerLen + kmerThreshold) for sparse envelope
#define MIN_KMERS_FOR_SPARSE_ENVELOPE 2

void DiagonalEnvelope::initFull() {
  LogThisAt(5, "Initializing full " << xLen << "*" << yLen << " envelope (no kmer-matching heuristic)" << endl);
  diagonals.clear();
  diagonals.reserve (xLen + yLen - 1);
  for (int d = minDiagonal(); d <= maxDiagonal(); ++d)
    diagonals.push_back (d);
  initStorage();
}

void DiagonalEnvelope::initSparse (const KmerIndex& yKmerIndex, unsigned int bandSize, int kmerThreshold, size_t cellSize, size_t maxSize) {
  const unsigned int kmerLen = yKmerIndex.kmerLen;
  
  if (kmerThreshold >= 0) {
    const SeqIdx minLenForSparse = MIN_KMERS_FOR_SPARSE_ENVELOPE * (kmerLen + kmerThreshold);
    if (px->length() < minLenForSparse || py->length() < minLenForSparse) {
      initFull();
      return;
    }
  }

  const vguard<AlphTok> xTok = px->tokens (yKmerIndex.alphabet);
  const AlphTok alphabetSize = (AlphTok) yKmerIndex.alphabet.size();
  
  map<int,unsigned int> diagKmerCount;
  for (SeqIdx i = 0; i <= xLen - kmerLen; ++i) {
    const auto yKmerIndexIter = yKmerIndex.kmerLocations.find (makeKmer (kmerLen, xTok.begin() + i, alphabetSize));
    if (yKmerIndexIter != yKmerIndex.kmerLocations.end())
      for (auto j : yKmerIndexIter->second)
	++diagKmerCount[get_diag(i,j)];
  }

  map<unsigned int,set<unsigned int> > countDistrib;
  for (const auto& diagKmerCountElt : diagKmerCount)
    countDistrib[diagKmerCountElt.second].insert (diagKmerCountElt.first);

  if (LoggingThisAt(7)) {
    LogStream (7, "Distribution of " << kmerLen << "-mer matches per diagonal for " << px->name << " vs " << py->name << ':' << endl);
    for (const auto& countDistribElt : countDistrib)
      LogStream (7, plural(countDistribElt.second.size(),"diagonal") << " have " << plural(countDistribElt.first,"match","matches") << endl);
  }

  set<int> diags, storageDiags;
  diags.insert(0);  // always add the zeroth diagonal to ensure at least one path exists
  storageDiags.insert(0);

  const unsigned int halfBandSize = bandSize / 2;
  const size_t diagSize = min(xLen,yLen) * cellSize;
  unsigned int nPastThreshold = 0;

  unsigned int threshold = numeric_limits<unsigned int>::max();
  bool foundThreshold = false;
  if (kmerThreshold >= 0) {
    threshold = kmerThreshold;
    foundThreshold = true;
  } else
    LogThisAt (5, "Automatically setting threshold based on memory limit of " << maxSize << " bytes (each diagonal takes " << diagSize << " bytes)" << endl);

  for (auto countDistribIter = countDistrib.crbegin();
       countDistribIter != countDistrib.crend();
       ++countDistribIter) {

    if (kmerThreshold >= 0 && countDistribIter->first < kmerThreshold)
      break;

    set<int> moreDiags = diags, moreStorageDiags = storageDiags;
    unsigned int moreNPastThreshold = nPastThreshold;
    for (auto seedDiag : countDistribIter->second) {
      ++moreNPastThreshold;
      const int dMin = max (minDiagonal(), (int) seedDiag - (int) halfBandSize);
      const int dMax = min (maxDiagonal(), (int) seedDiag + (int) halfBandSize);
      for (int d = dMin; d <= dMax; ++d)
	moreDiags.insert (d);
      for (int d = dMin - 1; d <= dMax + 1; ++d)
	moreStorageDiags.insert (d);
    }

    if (kmerThreshold < 0) {
      if (moreStorageDiags.size() * diagSize >= maxSize)
	break;
      threshold = countDistribIter->first;
      foundThreshold = true;
    }
    swap (diags, moreDiags);
    swap (storageDiags, moreStorageDiags);
    nPastThreshold = moreNPastThreshold;
  }

  if (foundThreshold)
    LogThisAt (5, "Threshold # of " << kmerLen << "-mer matches for seeding a diagonal is " << threshold << "; " << plural((long) nPastThreshold,"diagonal") << " over this threshold" << endl);
  else
    LogThisAt (5, "Couldn't find a suitable threshold that would fit within memory limit" << endl);
  LogThisAt (5, plural((long) diags.size(),"diagonal") << " in envelope (band size " << bandSize << "); estimated memory <" << (((storageDiags.size() * diagSize) >> 20) + 1) << "MB" << endl);

  diagonals = vguard<int> (diags.begin(), diags.end());
  initStorage();
}

void DiagonalEnvelope::initStorage() {
  set<int> storageDiags;
  for (auto d : diagonals) {
    storageDiags.insert (d);
    storageDiags.insert (d - 1);
    storageDiags.insert (d + 1);
  }
  storageDiagonals = vguard<int> (storageDiags.begin(), storageDiags.end());
  storageIndex = vguard<int> (xLen + yLen + 1, -1);
  for (size_t n = 0; n < storageDiagonals.size(); ++n)
    storageIndex[yLen + storageDiagonals[n]] = (int) n;
  storageOffset = vguard<int> (yLen + 1, -1);
  storageSize = vguard<size_t> (yLen + 1);
  cumulStorageSize = vguard<size_t> (yLen + 1);
  totalStorageSize = 0;
  for (SeqIdx j = 0; j <= yLen; ++j) {
    const vguard<int>::const_iterator b = storageBeginIntersecting(j);
    const vguard<int>::const_iterator e = storageEndIntersecting(j);
    storageSize[j] = e - b;
    cumulStorageSize[j] = totalStorageSize;
    totalStorageSize += storageSize[j];
    if (b != e)
      storageOffset[j] = storageIndex[yLen + *b];
  }
  LogThisAt(6, "Envelope for " << px->name << " vs " << py->name << " has " << totalStorageSize << " cells" << endl);
}

vguard<SeqIdx> DiagonalEnvelope::forward_i (SeqIdx j) const {
  vguard<SeqIdx> i_vec;
  i_vec.reserve (diagonals.size());
  for (auto d : diagonals)
    if (intersects (j, d))
      i_vec.push_back (get_i (j, d));
  return i_vec;
}

vguard<SeqIdx> DiagonalEnvelope::reverse_i (SeqIdx j) const {
  const vguard<SeqIdx> f = forward_i (j);
  return vguard<SeqIdx> (f.rbegin(), f.rend());
}
