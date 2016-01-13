#include <zlib.h>
#include <iostream>
#include "fastseq.h"
#include "util.h"
#include "logger.h"

KSEQ_INIT(gzFile, gzread)

const string dnaAlphabet ("ACGT");

UnvalidatedAlphTok tokenize (char c, const string& alphabet) {
  UnvalidatedAlphTok tok;
  const char* alphStr = alphabet.c_str(); 
  tok = (UnvalidatedAlphTok) (strchr (alphStr, toupper(c)) - alphStr);
  return tok >= (int) strlen(alphStr) ? -1 : tok;
}

AlphTok dnaComplement (AlphTok token) {
  return dnaAlphabetSize - 1 - token;
}

char dnaComplementChar (char c) {
  const int tok = tokenize (c, dnaAlphabet);
  return tok < 0 ? c : dnaAlphabet[dnaComplement(tok)];
}

Kmer makeKmer (SeqIdx k, vector<unsigned int>::const_iterator tok, AlphTok alphabetSize) {
  Kmer kmer = 0, mul = 1;
  for (SeqIdx j = 0; j < k; ++j) {
    const unsigned int token = tok[k - j - 1];
    kmer += mul * token;
    mul *= alphabetSize;
  }
  return kmer;
}

Kmer numberOfKmers (SeqIdx k, AlphTok alphabetSize) {
  Kmer n;
  for (n = 1; k > 0; --k)
    n *= alphabetSize;
  return n;
}

string kmerToString (Kmer kmer, SeqIdx k, const string& alphabet) {
  string rev;
  for (SeqIdx j = 0; j < k; ++j, kmer = kmer / alphabet.size())
    rev += alphabet[kmer % alphabet.size()];
  return string (rev.rbegin(), rev.rend());
}

SeqIntervalCoords SeqIntervalCoords::compose (const SeqIntervalCoords& srcCoords) const {
  if (srcCoords.isNull())
    return *this;
  SeqIntervalCoords coords;
  coords.name = srcCoords.name;
  coords.rev = rev != srcCoords.rev;
  if (srcCoords.rev) {
    coords.start = srcCoords.end - end + 1;
    coords.end = srcCoords.end - start + 1;
  } else {
    coords.start = start + srcCoords.start - 1;
    coords.end = end + srcCoords.start - 1;
  }
  return coords;
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
      cerr << "Unknown symbol " << c << " in sequence " << name << endl;
      throw;
    }
    tok.push_back (t);
  }
  return tok;
}

vguard<Kmer> FastSeq::kmers (const string& alphabet, unsigned int k) const {
  if (k == 0)
    return vguard<Kmer> (length(), 0);
  vguard<AlphTok> tok = tokens(alphabet);
  vguard<int> count (alphabet.size(), 0);
  for (auto t : tok)
    ++count[t];
  const AlphTok mostFrequentToken = (AlphTok) (max_element(count.begin(),count.end()) - count.begin());
  tok.insert (tok.begin(), k - 1, mostFrequentToken);
  vector<Kmer> result;
  result.reserve (length());
  for (SeqIdx pos = 0; pos < length(); ++pos)
    result.push_back (makeKmer (k, tok.begin() + pos, (AlphTok) alphabet.size()));
  return result;
}

vguard<QualScore> FastSeq::qualScores() const {
  vguard<QualScore> q;
  if (hasQual()) {
    q.reserve (length());
    for (const auto& c : qual)
      q.push_back (qualScoreForChar (c));
  }
  return q;
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
    const z_off_t filepos = gztell (fp);
    if (kseq_read(ks) == -1)
      break;

    FastSeq seq;
    initFastSeq (seq, ks);
    seq.filename = filename;
    seq.filepos = filepos;

    seqs.push_back (seq);
  }
  kseq_destroy (ks);
  gzclose (fp);

  LogThisAt(3, "Read " << plural(seqs.size(),"sequence") << " from " << filename << endl);
  
  if (seqs.empty())
    Warn ("Couldn't read any sequences from %s", filename);
  
  return seqs;
}

FastSeq readIndexedFastSeq (const char* filename, z_off_t filepos) {
  gzFile fp = gzopen(filename, "r");
  Require (fp != Z_NULL, "Couldn't open %s", filename);

  Require (gzseek (fp, filepos, SEEK_SET) == filepos, "Couldn't seek to byte %ld in %s", (long) filepos, filename);

  kseq_t *ks = kseq_init(fp);
  Require (kseq_read(ks) != -1, "Couldn't read sequence starting at byte %ld in %s", (long) filepos, filename);

  FastSeq seq;
  initFastSeq (seq, ks);
  seq.filename = filename;
  seq.filepos = filepos;

  kseq_destroy (ks);
  gzclose (fp);

  LogThisAt(3, "Read 1 sequence from " << filename << " at offset " << filepos << endl);

  return seq;
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

string revcomp (const string& dnaSeq) {
  string rev = dnaSeq;
  const size_t len = dnaSeq.size();
  for (size_t i = 0; i < len; ++i)
    rev[len - 1 - i] = dnaComplementChar (dnaSeq[i]);
  return rev;
}

FastSeq FastSeq::revcomp() const {
  FastSeq fs;
  fs.name = "revcomp(" + name + ")";
  fs.comment = comment;
  fs.seq = ::revcomp(seq);
  fs.qual = string (qual.rbegin(), qual.rend());
  fs.source.name = name;
  fs.source.start = 1;
  fs.source.end = length();
  fs.source.rev = true;
  fs.source = fs.source.compose (source);
  return fs;
}

void addRevcomps (vguard<FastSeq>& db) {
  vguard<FastSeq> revcomps;
  revcomps.reserve (db.size());
  for (const auto& fs : db)
    revcomps.push_back (fs.revcomp());
  db.insert (db.end(), revcomps.begin(), revcomps.end());
}

KmerIndex::KmerIndex (const FastSeq& seq, const string& alphabet, SeqIdx kmerLen)
  : seq(seq), alphabet(alphabet), kmerLen(kmerLen)
{
  LogThisAt(5, "Building " << kmerLen << "-mer index for " << seq.name << endl);
  const vguard<AlphTok> tok = seq.tokens (alphabet);
  const AlphTok alphabetSize = (AlphTok) alphabet.size();
  const SeqIdx seqLen = seq.length();
  for (SeqIdx j = 0; j <= seqLen - kmerLen; ++j)
    kmerLocations[makeKmer (kmerLen, tok.begin() + j, alphabetSize)].push_back (j);

  if (LoggingThisAt(8)) {
    LogStream (8, "Frequencies of " << kmerLen << "-mers in " << seq.name << ':' << endl);
    for (const auto& kl : kmerLocations) {
      LogStream (8, kmerToString (kl.first, kmerLen, alphabet) << ' ' << kl.second.size() << endl);
    }
  }
}

map<string,size_t> FastSeq::makeNameIndex (const vguard<FastSeq>& y) {
  map<string,size_t> yDict;
  for (size_t ny = 0; ny < y.size(); ++ny)
    yDict[y[ny].name] = ny;
  return yDict;
}
