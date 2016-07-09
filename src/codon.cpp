#include "codon.h"

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
}
