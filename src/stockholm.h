#ifndef STOCKHOLM_INCLUDED
#define STOCKHOLM_INCLUDED

#include <map>
#include "fastseq.h"
#include "tree.h"
#include "alignpath.h"

#define StockholmNewHampshireTag "NH"
#define StockholmIDTag "ID"
#define StockholmLogProbTag "LP"

#define MinStockholmCharsPerRow 10
#define DefaultStockholmRowLength 80

struct Stockholm {
  vguard<FastSeq> gapped;
  map<string,string> gc;  // gc[tag][col]
  map<string,vguard<string> > gf;  // gf[tag][line]
  map<string,map<string,string> > gr;  // gr[tag][seqname][col]
  map<string,map<string,vguard<string> > > gs;  // gs[tag][seqname][line]

  Stockholm();
  Stockholm (istream& in);
  Stockholm (const vguard<FastSeq>& seq);
  Stockholm (const vguard<FastSeq>& seq, const Tree& tree);
  
  void read (istream& in);
  void write (ostream& out, size_t charsPerRow = DefaultStockholmRowLength) const;  // set charsPerRow to zero to unlimit

  void setTree (const Tree& tree, const char* tag = StockholmNewHampshireTag);
  Tree getTree() const;
  bool hasTree() const;
  
  size_t rows() const;
  size_t columns() const;
  AlignPath path() const;

  void assertFlush() const { (void) columns(); }
};

#endif /* STOCKHOLM_INCLUDED */
