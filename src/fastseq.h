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

UnvalidatedAlphTok tokenize (char c, const string& alphabet);

struct FastSeq {
  // basic FASTQ data
  string name, comment, seq, qual;
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
  void writeFasta (ostream& out) const;
  void writeFastq (ostream& out) const;
};

vguard<FastSeq> readFastSeqs (const char* filename);
void writeFastaSeqs (ostream& out, const vguard<FastSeq>& fastSeqs);
void writeFastqSeqs (ostream& out, const vguard<FastSeq>& fastSeqs);

set<string> fastSeqDuplicateNames (const vguard<FastSeq>& seqs);

#endif /* KSEQCONTAINER_INCLUDED */
