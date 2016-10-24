#include "ctok.h"
#include "util.h"
#include "logger.h"

UniversalCodonTokenizer codonTokenizer;

void CodonTokenizer::addToken (char tok, const char* cod, bool isStop) {
  string cods (tolower (cod));
  tok2cod[tok] = cods;
  cod2tok[cods] = tok;
  for (auto& c : cods)
    if (c == 't')
      c = 'u';
  cod2tok[cods] = tok;
  if (isStop)
    stopTok.insert (tok);
}

void CodonTokenizer::addWildGap() {
  addToken('-',"---");
  addToken('*',"***");
}

UniversalCodonTokenizer::UniversalCodonTokenizer() {
  addToken('K',"aaa");
  addToken('n',"aac");
  addToken('k',"aag");
  addToken('N',"aat");
  addToken('~',"aca");
  addToken('t',"acc");
  addToken('`',"acg");
  addToken('T',"act");
  addToken('3',"aga");
  addToken('#',"agc");
  addToken(']',"agg");
  addToken('%',"agt");
  addToken('|',"ata");
  addToken('i',"atc");
  addToken('M',"atg");
  addToken('I',"att");
  addToken('Q',"caa");
  addToken('h',"cac");
  addToken('q',"cag");
  addToken('H',"cat");
  addToken(',',"cca");
  addToken('p',"ccc");
  addToken('8',"ccg");
  addToken('P',"cct");
  addToken('=',"cga");
  addToken('r',"cgc");
  addToken('}',"cgg");
  addToken('R',"cgt");
  addToken('{',"cta");
  addToken('[',"ctc");
  addToken('/',"ctg");
  addToken('<',"ctt");
  addToken('E',"gaa");
  addToken('d',"gac");
  addToken('e',"gag");
  addToken('D',"gat");
  addToken('4',"gca");
  addToken('a',"gcc");
  addToken('&',"gcg");
  addToken('A',"gct");
  addToken('9',"gga");
  addToken('g',"ggc");
  addToken('6',"ggg");
  addToken('G',"ggt");
  addToken('^',"gta");
  addToken('v',"gtc");
  addToken('7',"gtg");
  addToken('V',"gtt");
  addToken('0',"taa",true);
  addToken('y',"tac");
  addToken('1',"tag",true);
  addToken('Y',"tat");
  addToken('5',"tca");
  addToken('s',"tcc");
  addToken('$',"tcg");
  addToken('S',"tct");
  addToken('2',"tga",true);
  addToken('c',"tgc");
  addToken('W',"tgg");
  addToken('C',"tgt");
  addToken('L',"tta");
  addToken('f',"ttc");
  addToken('l',"ttg");
  addToken('F',"ttt");

  addWildGap();
}

string CodonTokenizer::tokenize (const string& gappedSeq, bool allowStopCodons, const char* name) const {
  Assert (gappedSeq.size() % 3 == 0, "Can't codon-tokenize %s: length (%d) is not a multiple of 3", name, gappedSeq.size());
  string tokSeq;
  tokSeq.reserve (gappedSeq.size() / 3);
  for (SeqIdx pos = 0; pos < gappedSeq.size(); pos += 3) {
    const string cod = tolower (gappedSeq.substr (pos, 3));
    Assert (cod2tok.count(cod), "Unknown codon '%s' at position %u in %s", cod.c_str(), pos, name);
    const char tok = cod2tok.at(cod);
    if (!allowStopCodons && isStopCodon(tok)) {
      if (pos + 3 == gappedSeq.size()) {
	LogThisAt (5, "Ignoring stop codon '" << cod << "' at end of " << name << endl);
	continue;
      }
      Abort ("Illegal stop codon '%s' at position %u in %s", cod.c_str(), pos, name);
    } else
      tokSeq.push_back (tok);
  }
  LogThisAt(7,"Tokenized " << name << " " << gappedSeq << " to " << tokSeq << endl);
  return tokSeq;
}

string CodonTokenizer::detokenize (const string& tokSeq) const {
  string gappedSeq;
  gappedSeq.reserve (gappedSeq.size() * 3);
  for (SeqIdx pos = 0; pos < tokSeq.size(); ++pos) {
    const char tok = tokSeq[pos];
    Assert (tok2cod.count(tok), "Can't detokenize '%c'", tok);
    gappedSeq.append (tok2cod.at(tok));
  }
  LogThisAt(7,"Detokenized " << tokSeq << " to " << gappedSeq << endl);
  return gappedSeq;
}

vguard<FastSeq> CodonTokenizer::tokenize (const vguard<FastSeq>& gappedSeq, bool allowStopCodons) const {
  vguard<FastSeq> tokSeq;
  tokSeq.reserve (gappedSeq.size());
  for (auto& fs: gappedSeq) {
    LogThisAt(6,"Tokenizing " << fs.name << endl);
    FastSeq tfs;
    tfs.name = fs.name;
    tfs.comment = fs.comment;
    tfs.seq = tokenize (fs.seq, allowStopCodons, fs.name.c_str());
    tokSeq.push_back (tfs);
  }
  return tokSeq;
}

Stockholm CodonTokenizer::tokenize (const Stockholm& stock, bool allowStopCodons) const {
  Stockholm tokStock;
  tokStock.gapped = tokenize (stock.gapped, allowStopCodons);
  tokStock.gf = stock.gf;
  tokStock.gs = stock.gs;
  return tokStock;
}

vguard<FastSeq> CodonTokenizer::detokenize (const vguard<FastSeq>& tokSeq) const {
  vguard<FastSeq> gappedSeq;
  gappedSeq.reserve (tokSeq.size());
  for (auto& tfs: tokSeq) {
    LogThisAt(6,"Untokenizing " << tfs.name << endl);
    FastSeq fs;
    fs.name = tfs.name;
    fs.comment = tfs.comment;
    fs.seq = detokenize (tfs.seq);
    gappedSeq.push_back (fs);
  }
  return gappedSeq;
}

Stockholm CodonTokenizer::detokenize (const Stockholm& tokStock) const {
  Stockholm stock;
  stock.gapped = detokenize (tokStock.gapped);
  stock.gf = tokStock.gf;
  stock.gs = tokStock.gs;
  return stock;
}

void CodonTokenizer::assertAlphabetTokenized (const string& alphabet) const {
  for (auto c: alphabet)
    Require (tok2cod.count(c), "Character %c in alphabet %s is not a tokenized codon", c, alphabet.c_str());
  if (alphabet.size() < 61)
    Warn ("Alphabet %s contains only %u characters, doesn't look like a tokenized codon alphabet", alphabet.c_str(), alphabet.size());
}

string CodonTokenizer::tokenAlphabet (bool allowStopCodons) const {
  const string dna ("tcag");
  string cod ("xxx"), alph;
  alph.reserve (allowStopCodons ? 64 : 61);
  for (size_t i = 0; i < 4; ++i) {
    cod[0] = dna[i];
    for (size_t j = 0; j < 4; ++j) {
      cod[1] = dna[j];
      for (size_t k = 0; k < 4; ++k) {
	cod[2] = dna[k];
	const char tok = cod2tok.at(cod);
	if (allowStopCodons || !isStopCodon(tok))
	  alph.push_back (tok);
      }
    }
  }
  return alph;
}
