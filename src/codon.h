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

  string tokenize (const string& gappedSeq) const;
  vguard<FastSeq> tokenize (const vguard<FastSeq>& gappedSeq) const;
  Stockholm tokenize (const Stockholm& stock) const;

  string untokenize (const string& gappedSeq) const;
  vguard<FastSeq> untokenize (const vguard<FastSeq>& gappedSeq) const;
  Stockholm untokenize (const Stockholm& stock) const;
};

extern CodonTokenizer codonTokenizer;

#endif /* CODON_INCLUDED */
