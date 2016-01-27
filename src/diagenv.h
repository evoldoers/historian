#ifndef DIAG_ENV_INCLUDED
#define DIAG_ENV_INCLUDED

#include <map>
#include "fastseq.h"
#include "util.h"

using namespace std;

// Somewhat arbitrary numbers:
// Default k-mer length
#define DEFAULT_KMER_LENGTH 6

// The number of kmer matches above which a diagonal is selected for a band
#define DEFAULT_KMER_THRESHOLD 14

// Default band size
#define DEFAULT_BAND_SIZE 64

// uncomment to add extra guards to storage indexing
// #define USE_DIAGONAL_ENVELOPE_GUARDS

struct DiagonalEnvelope {
  const FastSeq *px, *py;
  const SeqIdx xLen, yLen;
  vguard<int> diagonals, storageDiagonals;   // Sorted ascending. (i,j) is on diagonal d if i-j=d
  vguard<int> storageIndex;  // storageIndex[yLen + storageDiagonals[n]] = n, or -1 if outside envelope
  vguard<int> storageOffset;  // storageOffset[j] = storageIndex of first diagonal intersecting j
  vguard<size_t> storageSize;  // storageSize[j] = number of diagonals intersecting j
  vguard<size_t> cumulStorageSize;  // cumulStorageSize[j] = sum_{j'=0}^{j-1} storageSize[j]
  size_t totalStorageSize;  // sum_{j=0}^{yLen} storageSize[j]
  DiagonalEnvelope (const FastSeq& x, const FastSeq& y)
    : px(&x), py(&y), xLen(px->length()), yLen(py->length()) { }
  void initFull();
  void initSparse (const KmerIndex& yKmerIndex,
		   unsigned int bandSize = DEFAULT_BAND_SIZE,
		   int kmerThreshold = DEFAULT_KMER_THRESHOLD,  // negative => use memory guides
		   size_t cellSize = sizeof(double),
		   size_t maxSize = 0);
  void initStorage();
  inline int getStorageIndexSafe (SeqIdx i, SeqIdx j) const {
    const int idx = storageIndex[yLen + i - j];
    const int offsetIdx = idx - storageOffset[j];
    return (idx < 0 || offsetIdx < 0 || offsetIdx >= (int) storageSize[j]) ? -1 : (offsetIdx + (int) cumulStorageSize[j]);
  }
  inline int getStorageIndexUnsafe (SeqIdx i, SeqIdx j) const {
#ifdef USE_DIAGONAL_ENVELOPE_GUARDS
    const int idx = getStorageIndexSafe (i, j);
    if (idx < 0)
      Abort ("Access error at (i=%u,j=%u) storageIndex=%d storageOffset=%d storageSize=%lu", i, j, storageIndex[yLen+i-j], storageOffset[j], storageSize[j]);
#else /* USE_DIAGONAL_ENVELOPE_GUARDS */
    const int idx = storageIndex[(int) yLen + i - j] - storageOffset[j] + (int) cumulStorageSize[j];
#endif /* USE_DIAGONAL_ENVELOPE_GUARDS */
    return idx;
  }
  inline static int getDiagonal (SeqIdx i, SeqIdx j) { return i - j; }
  inline int minDiagonal() const { return 1 - (int) yLen; }
  inline int maxDiagonal() const { return xLen - 1; }
  inline int minStorageDiagonal() const { return - (int) yLen; }
  inline int maxStorageDiagonal() const { return xLen; }
  inline bool contains (SeqIdx i, SeqIdx j) const {
    const int diag = i - j;
    const vguard<int>::const_iterator iter = lower_bound (diagonals.begin(), diagonals.end(), diag);
    return iter != diagonals.end() && *iter == diag;
  }
  inline bool intersects (SeqIdx j, int diag) const {
    const int i = diag + j;
    return i > 0 && i <= xLen;
  }
  static inline SeqIdx get_i (SeqIdx j, int diag) { return (SeqIdx) (diag + j); }
  static inline int get_diag (SeqIdx i, SeqIdx j) { return (int) i - (int) j; }
  vguard<SeqIdx> forward_i (SeqIdx j) const;
  vguard<SeqIdx> reverse_i (SeqIdx j) const;

  inline vguard<int>::const_iterator beginIntersecting (SeqIdx j) const {
    // i = diag + j
    // want i > 0 && i <= xLen
    // so diag > -j && diag <= xLen - j
    return upper_bound (diagonals.begin(), diagonals.end(), -(int)j);
  }

  inline vguard<int>::const_iterator endIntersecting (SeqIdx j) const {
    // see comment in beginIntersecting()
    return upper_bound (diagonals.begin(), diagonals.end(), (int) xLen - (int) j);
  }

  inline vguard<int>::const_iterator storageBeginIntersecting (SeqIdx j) const {
    // i = diag + j
    // want i >= 0 (for storage, as opposed to iterator)
    // so diag > -j-1 && diag <= xLen - j
    return upper_bound (storageDiagonals.begin(), storageDiagonals.end(), -1 - (int)j);
  }

  inline vguard<int>::const_iterator storageEndIntersecting (SeqIdx j) const {
    return upper_bound (storageDiagonals.begin(), storageDiagonals.end(), (int) xLen - (int) j);
  }

  class DiagonalEnvelopeIterator : public iterator<forward_iterator_tag,SeqIdx> {
  private:
    const SeqIdx j;
    vguard<int>::const_iterator iter;
    const vguard<int>::const_iterator iterEnd;
  public:
    inline DiagonalEnvelopeIterator (SeqIdx j, const vguard<int>::const_iterator& iter, const vguard<int>::const_iterator& iterEnd)
      : j(j), iter(iter), iterEnd(iterEnd)
    { }
    inline bool operator== (const DiagonalEnvelopeIterator& dei) const { return iter == dei.iter; }
    inline bool operator!= (const DiagonalEnvelopeIterator& dei) const { return iter != dei.iter; }
    inline SeqIdx operator*() const { return (SeqIdx) (*iter + j); }
    inline DiagonalEnvelopeIterator& operator++() { ++iter; return *this; }
    inline DiagonalEnvelopeIterator operator++(int) { DiagonalEnvelopeIterator dei = *this; ++*this; return dei; }
    inline bool finished() const { return iter == iterEnd; }
  };

  class DiagonalEnvelopeReverseIterator : public iterator<forward_iterator_tag,SeqIdx> {
  private:
    const SeqIdx j;
    vguard<int>::const_iterator iter;
    const vguard<int>::const_iterator iterEnd;
  public:
    inline DiagonalEnvelopeReverseIterator (SeqIdx j, const vguard<int>::const_iterator& iter, const vguard<int>::const_iterator& iterEnd)
      : j(j), iter(iter), iterEnd(iterEnd)
    { }
    inline bool operator== (const DiagonalEnvelopeReverseIterator& dei) const { return iter == dei.iter; }
    inline bool operator!= (const DiagonalEnvelopeReverseIterator& dei) const { return iter != dei.iter; }
    inline SeqIdx operator*() const { return (SeqIdx) (*(iter - 1) + j); }
    inline DiagonalEnvelopeReverseIterator& operator++() { --iter; return *this; }
    inline DiagonalEnvelopeReverseIterator operator++(int) { DiagonalEnvelopeReverseIterator deri = *this; ++*this; return deri; }
    inline bool finished() const { return iter == iterEnd; }
  };

  typedef DiagonalEnvelopeIterator iterator;
  typedef const DiagonalEnvelopeIterator const_iterator;
  typedef DiagonalEnvelopeReverseIterator reverse_iterator;
  typedef const DiagonalEnvelopeReverseIterator const_reverse_iterator;

  inline iterator begin (SeqIdx j) const { return DiagonalEnvelopeIterator (j, beginIntersecting(j), endIntersecting(j)); }
  inline const_iterator end (SeqIdx j) const { return DiagonalEnvelopeIterator (j, endIntersecting(j), endIntersecting(j)); }

  inline reverse_iterator rbegin (SeqIdx j) const { return DiagonalEnvelopeReverseIterator (j, endIntersecting(j), beginIntersecting(j)); }
  inline const_reverse_iterator rend (SeqIdx j) const { return DiagonalEnvelopeReverseIterator (j, beginIntersecting(j), beginIntersecting(j)); }
};

#endif /* DIAG_ENV_INCLUDED */
