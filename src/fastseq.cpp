#include <zlib.h>
#include <iostream>
#include "fastseq.h"
#include "util.h"
#include "logger.h"

KSEQ_INIT(gzFile, gzread)

UnvalidatedAlphTok tokenize (char c, const string& alphabet) {
  const char* alphStr = alphabet.c_str();
  const char* ptok = strchr (alphStr, c);
  if (ptok == NULL)
    ptok = strchr (alphStr, isupper(c) ? tolower(c) : toupper(c));
  return ptok ? (UnvalidatedAlphTok) (ptok - alphStr) : -1;
}

const char FastSeq::minQualityChar = '!';
const char FastSeq::maxQualityChar = '~';
const QualScore FastSeq::qualScoreRange = 94;

vguard<AlphTok> FastSeq::tokens (const string& alphabet) const {
  vguard<AlphTok> tok;
  tok.reserve (length());
  for (const auto& c : seq) {
    const int t = tokenize (c, alphabet);
    if (t < 0) {
      cerr << "Unknown symbol " << c << " in sequence " << name << " (alphabet is " << alphabet << ")" << endl;
      throw;
    }
    tok.push_back (t);
  }
  return tok;
}

void FastSeq::writeFasta (ostream& out) const {
  out << '>' << name;
  if (comment.size())
    out << ' ' << comment;
  out << endl;
  out << seq << endl;
}

void FastSeq::writeFastq (ostream& out) const {
  out << '@' << name;
  if (comment.size())
    out << ' ' << comment;
  out << endl;
  out << seq << endl;
  if (hasQual())
    out << '+' << endl << qual << endl;
}

void writeFastaSeqs (ostream& out, const vguard<FastSeq>& fastSeqs) {
  for (const auto& s : fastSeqs)
    s.writeFasta (out);
}

void writeFastqSeqs (ostream& out, const vguard<FastSeq>& fastSeqs) {
  for (const auto& s : fastSeqs)
    s.writeFastq (out);
}

void initFastSeq (FastSeq& seq, kseq_t* ks) {
  seq.name = string(ks->name.s);
  seq.seq = string(ks->seq.s);
  if (ks->comment.l)
    seq.comment = string(ks->comment.s);
  if (ks->qual.l == ks->seq.l)
    seq.qual = string(ks->qual.s);
}

vguard<FastSeq> readFastSeqs (const char* filename) {
  vguard<FastSeq> seqs;

  gzFile fp = gzopen(filename, "r");
  Require (fp != Z_NULL, "Couldn't open %s", filename);

  kseq_t *ks = kseq_init(fp);
  while (true) {
    if (kseq_read(ks) == -1)
      break;

    FastSeq seq;
    initFastSeq (seq, ks);

    seqs.push_back (seq);
  }
  kseq_destroy (ks);
  gzclose (fp);

  LogThisAt(3, "Read " << plural(seqs.size(),"sequence") << " from " << filename << endl);
  
  if (seqs.empty())
    Warn ("Couldn't read any sequences from %s", filename);
  
  return seqs;
}

set<string> fastSeqDuplicateNames (const vguard<FastSeq>& seqs) {
  set<string> name, dups;
  for (const auto& s : seqs) {
    if (name.find(s.name) != name.end())
      dups.insert (s.name);
    name.insert (s.name);
  }
  return dups;
}
