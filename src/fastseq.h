#ifndef FASTSEQ_INCLUDED
#define FASTSEQ_INCLUDED

#include <string>
#include <algorithm>
#include <map>
#include <set>
#include <zlib.h>
#include "vguard.h"
#include "kseq.h"

using namespace std;

#define DefaultFastaCharsPerLine 50
#define InvalidAlphabetToken -1

// alphabets
typedef unsigned int SeqIdx;
typedef unsigned int AlphTok;
typedef int UnvalidatedAlphTok;
typedef unsigned long long Kmer;
typedef unsigned int QualScore;
typedef string AlphSym;

typedef vguard<AlphTok> TokSeq;
typedef vguard<UnvalidatedAlphTok> UnvalidatedTokSeq;

UnvalidatedAlphTok tokenize (char c, const string& alphabet);  // returns -1 if invalid
UnvalidatedTokSeq tokenize (const string& s, const string& alphabet);
TokSeq validTokenize (const string& s, const string& alphabet, const char* seqname = NULL);
string detokenize (const TokSeq& s, const string& alphabet);

// kmers
bool kmerValid (SeqIdx k, vector<UnvalidatedAlphTok>::const_iterator tok);
Kmer makeKmer (SeqIdx k, vector<UnvalidatedAlphTok>::const_iterator tok, AlphTok alphabetSize);
Kmer numberOfKmers (SeqIdx k, AlphTok alphabetSize);
string kmerToString (Kmer kmer, SeqIdx k, const string& alphabet);

// extended alphabets
typedef map<AlphSym,AlphTok> ExtendedAlphabet;
typedef map<AlphTok,AlphSym> ExtendedAlphabetIndex;
size_t minSymbolLength (const ExtendedAlphabet& alphabet);
ExtendedAlphabetIndex makeAlphabetIndex (const ExtendedAlphabet& alphabet);

UnvalidatedTokSeq tokenize (const string& s, const ExtendedAlphabet& alphabet);
TokSeq validTokenize (const string& s, const ExtendedAlphabet& alphabet, const char* seqname = NULL);
string detokenize (const TokSeq& s, const ExtendedAlphabet& alphabet);

// FASTA- and FASTQ-format sequences
struct FastSeq {
  // basic FASTQ data
  string name, comment, seq, qual;
  // statics
  static const char minQualityChar, maxQualityChar;
  static const QualScore qualScoreRange;
  // methods
  SeqIdx length() const { return (SeqIdx) seq.size(); }
  bool hasQual() const { return qual.size() == length(); }
  static inline QualScore qualScoreForChar (char c) {
    return max (0, min ((int) qualScoreRange - 1, (int) (c - minQualityChar)));
  }
  static char charForQualScore (QualScore q) {
    return max (minQualityChar, min (maxQualityChar, (char) (q + minQualityChar)));
  }
  inline QualScore getQualScoreAt (SeqIdx pos) const { return qualScoreForChar (qual[pos]); }
  TokSeq tokens (const string& alphabet) const;
  UnvalidatedTokSeq unvalidatedTokens (const string& alphabet) const;
  void writeFasta (ostream& out) const;
  void writeFastq (ostream& out) const;
};

vguard<FastSeq> readFastSeqs (const char* filename);
void writeFastaSeqs (ostream& out, const vguard<FastSeq>& fastSeqs);
void writeFastqSeqs (ostream& out, const vguard<FastSeq>& fastSeqs);

set<string> fastSeqDuplicateNames (const vguard<FastSeq>& seqs);

struct KmerIndex {
  const FastSeq& seq;
  const string& alphabet;
  const SeqIdx kmerLen;
  map <Kmer, vector<SeqIdx> > kmerLocations;
  KmerIndex (const FastSeq& seq, const string& alphabet, SeqIdx kmerLen);
};

#endif /* KSEQCONTAINER_INCLUDED */
