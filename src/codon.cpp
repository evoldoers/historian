#include "codon.h"
#include "util.h"

CodonTokenizer codonTokenizer;

void CodonTokenizer::addToken (char tok, const char* cod) {
  const string cods (cod);
  tok2cod[tok] = cods;
  cod2tok[cods] = tok;
}

CodonTokenizer::CodonTokenizer() {
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
  addToken('0',"taa"); 
  addToken('y',"tac"); 
  addToken('1',"tag"); 
  addToken('Y',"tat"); 
  addToken('5',"tca"); 
  addToken('s',"tcc"); 
  addToken('$',"tcg"); 
  addToken('S',"tct"); 
  addToken('2',"tga"); 
  addToken('c',"tgc"); 
  addToken('W',"tgg"); 
  addToken('C',"tgt"); 
  addToken('L',"tta"); 
  addToken('f',"ttc"); 
  addToken('l',"ttg"); 
  addToken('F',"ttt");   

  addToken('-',"---");
  addToken('*',"***");
}

string CodonTokenizer::tokenize (const string& gappedSeq) const {
  Assert (gappedSeq.size() % 3 == 0, "Can't codon-tokenize a sequence that's not a multiple of 3 nucleotides in length");
  string tokSeq;
  tokSeq.reserve (gappedSeq.size() / 3);
  for (SeqIdx pos = 0; pos < gappedSeq.size(); pos += 3) {
    const string cod = tolower (gappedSeq.substr (pos, 3));
    Assert (cod2tok.count(cod), "Can't tokenize '%s'", cod.c_str());
    tokSeq.push_back (cod2tok.at(cod));
  }
  return tokSeq;
}

string CodonTokenizer::untokenize (const string& tokSeq) const {
  string gappedSeq;
  gappedSeq.reserve (gappedSeq.size() * 3);
  for (SeqIdx pos = 0; pos < tokSeq.size(); ++pos) {
    const char tok = tokSeq[pos];
    Assert (tok2cod.count(tok), "Can't untokenize '%c'", tok);
    gappedSeq.append (tok2cod.at(tok));
  }
  return gappedSeq;
}

vguard<FastSeq> CodonTokenizer::tokenize (const vguard<FastSeq>& gappedSeq) const {
  vguard<FastSeq> tokSeq;
  tokSeq.reserve (gappedSeq.size());
  for (auto& fs: gappedSeq) {
    FastSeq tfs;
    tfs.name = fs.name;
    tfs.comment = fs.comment;
    tfs.seq = tokenize (fs.seq);
    tokSeq.push_back (tfs);
  }
  return tokSeq;
}

Stockholm CodonTokenizer::tokenize (const Stockholm& stock) const {
  Stockholm tokStock;
  tokStock.gapped = tokenize (stock.gapped);
  tokStock.gf = stock.gf;
  tokStock.gs = stock.gs;
  return tokStock;
}

vguard<FastSeq> CodonTokenizer::untokenize (const vguard<FastSeq>& tokSeq) const {
  vguard<FastSeq> gappedSeq;
  gappedSeq.reserve (tokSeq.size());
  for (auto& tfs: tokSeq) {
    FastSeq fs;
    fs.name = tfs.name;
    fs.comment = tfs.comment;
    fs.seq = untokenize (tfs.seq);
    gappedSeq.push_back (fs);
  }
  return gappedSeq;
}

Stockholm CodonTokenizer::untokenize (const Stockholm& tokStock) const {
  Stockholm stock;
  stock.gapped = untokenize (tokStock.gapped);
  stock.gf = tokStock.gf;
  stock.gs = tokStock.gs;
  return stock;
}
