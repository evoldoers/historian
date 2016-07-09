#ifndef CODON_INCLUDED
#define CODON_INCLUDED

#include "fastseq.h"
#include "stockholm.h"

struct CodonTokenizer
{
  map<string,char> cod2tok;
  map<char,string> tok2cod;

  CodonTokenizer();

  vguard<FastSeq> tokenize (const vguard<FastSeq>& gappedSeq) const;
  Stockholm tokenize (const Stockholm& stock) const;

  vguard<FastSeq> untokenize (const vguard<FastSeq>& gappedSeq) const;
  Stockholm untokenize (const Stockholm& stock) const;
};

extern CodonTokenizer codonTokenizer;

#endif /* CODON_INCLUDED */
