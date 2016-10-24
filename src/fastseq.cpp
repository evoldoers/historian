#include <zlib.h>
#include <iostream>
#include "fastseq.h"
#include "util.h"
#include "alignpath.h"
#include "logger.h"

KSEQ_INIT(gzFile, gzread)

UnvalidatedAlphTok tokenize (char c, const string& alphabet) {
  const char* alphStr = alphabet.c_str();
  const char* ptok = strchr (alphStr, c);
  if (ptok == NULL)
    ptok = strchr (alphStr, isupper(c) ? tolower(c) : toupper(c));
  return ptok ? (UnvalidatedAlphTok) (ptok - alphStr) : InvalidAlphabetToken;
}

UnvalidatedTokSeq tokenize (const string& s, const string& alphabet) {
  UnvalidatedTokSeq tok;
  tok.reserve (s.size());
  for (auto c: s) {
    const auto t = tokenize (c, alphabet);
    tok.push_back (t);
  }
  return tok;
}

const char FastSeq::minQualityChar = '!';
const char FastSeq::maxQualityChar = '~';
const QualScore FastSeq::qualScoreRange = 94;

TokSeq validTokenize (const string& s, const string& alphabet, const char* seqname) {
  TokSeq tok;
  tok.reserve (s.size());
  vguard<int> freq (alphabet.size());
  for (auto c: s) {
    const auto t = tokenize (c, alphabet);
    if (t >= 0)
      ++freq[t];
  }
  const AlphTok mostFrequentTok = (AlphTok) (max_element (freq.begin(), freq.end()) - freq.begin());
  for (auto c: s) {
    UnvalidatedAlphTok t = Alignment::isWildcard(c) ? mostFrequentTok : tokenize (c, alphabet);
    if (t < 0) {
      cerr << "Unknown symbol " << c << " in sequence"
	   << (seqname ? ((string(" ") + seqname)) : string())
	   << " (alphabet is " << alphabet << "). Replacing with most frequent symbol "
	   << alphabet[mostFrequentTok] << endl;
      t = mostFrequentTok;
      //      throw;
    }
    tok.push_back ((AlphTok) t);
  }
  return tok;
}

string detokenize (const TokSeq& toks, const string& alphabet) {
  string seq;
  seq.reserve (toks.size());
  for (auto tok : toks) {
    if (tok >= alphabet.size()) {
      cerr << "Unknown token " << tok << " in sequence (alphabet is " << alphabet << ")" << endl;
      throw;
    }
    seq.push_back (alphabet[tok]);
  }
  return seq;
}

size_t minSymbolLength (const ExtendedAlphabet& alphabet) {
  vguard<size_t> len;
  len.reserve (alphabet.size());
  for (const auto& sym: extract_keys(alphabet))
    len.push_back (sym.size());
  return *min_element (len.begin(), len.end());
}

ExtendedAlphabetIndex makeAlphabetIndex (const ExtendedAlphabet& alphabet) {
  ExtendedAlphabetIndex inv;
  for (const auto& sym_tok: alphabet)
    inv[sym_tok.second] = sym_tok.first;
  return inv;
}

UnvalidatedTokSeq tokenize (const string& s, const ExtendedAlphabet& alphabet) {
  UnvalidatedTokSeq result;
  const size_t minLen = minSymbolLength(alphabet);
  size_t pos = 0;
  while (pos < s.size()) {
    UnvalidatedAlphTok tok = InvalidAlphabetToken;
    size_t len;
    for (len = minLen; pos + len <= s.size(); ++len) {
      const AlphSym sym = s.substr (pos, len);
      if (alphabet.count(sym)) {
	tok = alphabet.at(sym);
	break;
      }
    }
    result.push_back (tok);
    pos += (tok < 0 ? 1 : len);
  }
  return result;
}

TokSeq validTokenize (const string& s, const ExtendedAlphabet& alphabet, const char* seqname) {
  TokSeq result;
  const size_t minLen = minSymbolLength(alphabet);
  size_t pos = 0;
  while (pos < s.size()) {
    bool found = false;
    size_t len;
    for (len = minLen; pos + len <= s.size(); ++len) {
      const AlphSym sym = s.substr (pos, len);
      if (alphabet.count(sym)) {
	result.push_back (alphabet.at(sym));
	found = true;
	break;
      }
    }
    if (!found)
      cerr << "Unknown token '" << s.substr(pos,minLen) << "' in sequence"
	   << (seqname ? ((string(" ") + seqname)) : string()) << endl;
    pos += found ? len : 1;
  }
  return result;
}

string detokenize (const TokSeq& s, const ExtendedAlphabet& alphabet) {
  string result;
  const size_t minLen = minSymbolLength(alphabet);
  result.reserve (minLen * s.size());
  const auto inv = makeAlphabetIndex (alphabet);
  for (auto tok: s)
    result.append (inv.count(tok) ? inv.at(tok) : string(1,Alignment::wildcardChar));
  return result;
}

TokSeq FastSeq::tokens (const string& alphabet) const {
  return validTokenize (seq, alphabet, name.c_str());
}

UnvalidatedTokSeq FastSeq::unvalidatedTokens (const string& alphabet) const {
  return tokenize (seq, alphabet);
}

bool kmerValid (SeqIdx k, vector<int>::const_iterator tok) {
  for (SeqIdx j = 0; j < k; ++j)
    if (tok[j] < 0)
      return false;
  return true;
}

Kmer makeKmer (SeqIdx k, vector<int>::const_iterator tok, AlphTok alphabetSize) {
  Kmer kmer = 0, mul = 1;
  for (SeqIdx j = 0; j < k; ++j) {
    const int token = tok[k - j - 1];
    Assert (token >= 0, "Invalid token in makeKmer");
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

void FastSeq::writeFasta (ostream& out) const {
  out << '>' << name;
  if (comment.size())
    out << ' ' << comment;
  out << endl;
  const size_t width = DefaultFastaCharsPerLine;
  for (size_t i = 0; i < seq.size(); i += width)
    out << seq.substr(i,width) << endl;
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
  if (ks->name.l)
    seq.name = string(ks->name.s);
  if (ks->seq.l)
    seq.seq = string(ks->seq.s);
  if (ks->comment.l)
    seq.comment = string(ks->comment.s);
  if (ks->qual.l && ks->qual.l == ks->seq.l)
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

KmerIndex::KmerIndex (const FastSeq& seq, const string& alphabet, SeqIdx kmerLen)
  : seq(seq), alphabet(alphabet), kmerLen(kmerLen)
{
  LogThisAt(5, "Building " << kmerLen << "-mer index for " << seq.name << endl);
  const UnvalidatedTokSeq tok = seq.unvalidatedTokens (alphabet);
  const AlphTok alphabetSize = (AlphTok) alphabet.size();
  const SeqIdx seqLen = seq.length();
  for (SeqIdx j = 0; j + kmerLen <= seqLen; ++j)
    if (kmerValid(kmerLen,tok.begin() + j))
      kmerLocations[makeKmer (kmerLen, tok.begin() + j, alphabetSize)].push_back (j);

  if (LoggingThisAt(8)) {
    LogStream (8, "Frequencies of " << kmerLen << "-mers in " << seq.name << ':' << endl);
    for (const auto& kl : kmerLocations) {
      LogStream (8, kmerToString (kl.first, kmerLen, alphabet) << ' ' << kl.second.size() << endl);
    }
  }
}
