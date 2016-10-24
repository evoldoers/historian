#ifndef CODON_INCLUDED
#define CODON_INCLUDED

#include "fastseq.h"
#include "stockholm.h"

class CodonTokenizer
{
private:
  map<string,char> cod2tok;
  map<char,string> tok2cod;

  set<char> stopTok;
  
protected:
  void addToken (char tok, const char* cod, bool isStop = false);
  void addWildGap();
  
public:
  char getToken (const string& codon) const { return cod2tok.at(codon); }
  const string& getCodon (const char tok) const { return tok2cod.at(tok); }
  
  string tokenize (const string& gappedSeq, bool allowStopCodons = false, const char* name = "sequence") const;
  vguard<FastSeq> tokenize (const vguard<FastSeq>& gappedSeq, bool allowStopCodons = false) const;
  Stockholm tokenize (const Stockholm& stock, bool allowStopCodons = false) const;

  string detokenize (const string& gappedSeq) const;
  vguard<FastSeq> detokenize (const vguard<FastSeq>& gappedSeq) const;
  Stockholm detokenize (const Stockholm& stock) const;

  string tokenAlphabet (bool allowStopCodons = false) const;
  void assertAlphabetTokenized (const string& alphabet) const;
  
  inline bool isStopCodon (const char c) const {
    return stopTok.count(c) > 0;
  };
};

class UniversalCodonTokenizer : public CodonTokenizer {
public:
  UniversalCodonTokenizer();
};

extern UniversalCodonTokenizer codonTokenizer;

#endif /* CODON_INCLUDED */
