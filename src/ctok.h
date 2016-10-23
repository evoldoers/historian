#ifndef CODON_INCLUDED
#define CODON_INCLUDED

#include "fastseq.h"
#include "stockholm.h"

class CodonTokenizer
{
private:
  map<string,char> cod2tok;
  map<char,string> tok2cod;

  void addToken (char tok, const char* cod);
  
public:
  CodonTokenizer();

  string tokenize (const string& gappedSeq, bool allowStopCodons = false, const char* name = "sequence") const;
  vguard<FastSeq> tokenize (const vguard<FastSeq>& gappedSeq, bool allowStopCodons = false) const;
  Stockholm tokenize (const Stockholm& stock, bool allowStopCodons = false) const;

  string detokenize (const string& gappedSeq) const;
  vguard<FastSeq> detokenize (const vguard<FastSeq>& gappedSeq) const;
  Stockholm detokenize (const Stockholm& stock) const;

  void assertAlphabetTokenized (const string& alphabet) const;
  
  inline static bool isStopCodon (const char c) {
    return c == '0' || c == '1' || c == '2';
  };
};

extern CodonTokenizer codonTokenizer;

#endif /* CODON_INCLUDED */
