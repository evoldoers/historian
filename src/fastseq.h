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

typedef unsigned int SeqIdx;
typedef unsigned int AlphTok;
typedef int UnvalidatedAlphTok;
typedef unsigned long long Kmer;
typedef unsigned int QualScore;

UnvalidatedAlphTok tokenize (char c, const string& alphabet);
Kmer makeKmer (SeqIdx k, vector<AlphTok>::const_iterator tok, AlphTok alphabetSize);
Kmer numberOfKmers (SeqIdx k, AlphTok alphabetSize);
string kmerToString (Kmer kmer, SeqIdx k, const string& alphabet);

#define dnaAlphabetSize 4
extern const string dnaAlphabet;
AlphTok dnaComplement (AlphTok token);
char dnaComplementChar (char c);

struct SeqIntervalCoords {
  string name;
  SeqIdx start, end;  // 1-based, closed (includes begin & end)
  bool rev;
  SeqIntervalCoords() : name(), start(0), end(0), rev(false) { }
  SeqIntervalCoords (const string& n, SeqIdx s, SeqIdx e, bool r)
    : name(n), start(s), end(e), rev(r)
  { }
  bool isNull() const { return name.empty(); }
  SeqIntervalCoords compose (const SeqIntervalCoords& srcCoords) const;
};

struct FastSeq {
  // basic FASTQ data
  string name, comment, seq, qual;
  // metadata
  SeqIntervalCoords source;  // if !source.isNull(), describes origin
  string filename;
  z_off_t filepos;
  // statics
  static const char minQualityChar, maxQualityChar;
  static const QualScore qualScoreRange;
  // methods
  FastSeq() : filepos(-1) { }
  SeqIdx length() const { return (SeqIdx) seq.size(); }
  bool hasQual() const { return qual.size() == length(); }
  static inline QualScore qualScoreForChar (char c) {
    return max (0, min ((int) qualScoreRange - 1, (int) (c - minQualityChar)));
  }
  static char charForQualScore (QualScore q) {
    return max (minQualityChar, min (maxQualityChar, (char) (q + minQualityChar)));
  }
  inline QualScore getQualScoreAt (SeqIdx pos) const { return qualScoreForChar (qual[pos]); }
  vguard<AlphTok> tokens (const string& alphabet) const;
  vguard<Kmer> kmers (const string& alphabet, unsigned int k) const;
  vguard<QualScore> qualScores() const;
  void writeFasta (ostream& out) const;
  void writeFastq (ostream& out) const;
  FastSeq revcomp() const;
  static map<string,size_t> makeNameIndex (const vguard<FastSeq>& y);
};

vguard<FastSeq> readFastSeqs (const char* filename);
FastSeq readIndexedFastSeq (const char* filename, z_off_t filepos);
void writeFastaSeqs (ostream& out, const vguard<FastSeq>& fastSeqs);
void writeFastqSeqs (ostream& out, const vguard<FastSeq>& fastSeqs);

set<string> fastSeqDuplicateNames (const vguard<FastSeq>& seqs);

string revcomp (const string& dnaSeq);
void addRevcomps (vguard<FastSeq>& db);

struct KmerIndex {
  const FastSeq& seq;
  const string& alphabet;
  const SeqIdx kmerLen;
  map <Kmer, vector<SeqIdx> > kmerLocations;
  KmerIndex (const FastSeq& seq, const string& alphabet, SeqIdx kmerLen);
};

#endif /* KSEQCONTAINER_INCLUDED */
