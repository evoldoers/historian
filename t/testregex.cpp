#include <iostream>
#include <fstream>
#include <string.h>
#include "../src/regexmacros.h"

using namespace std;

#define Q(x) #x
#define QUOTE(x) Q(x)
#define TEST2(X,Y) if (!(X)) { cout << "failed " << QUOTE(X) << Y << "\n"; exit (EXIT_FAILURE); }
#define TEST(X) TEST2(X,"")

bool test_regex (const regex& re, const char* c1, const char* c2) {
  string s1(c1), s2(c2);
  smatch sm;
  if (regex_match (s2, sm, re))
    return false;
  return regex_match (s1, sm, re);
}

int main (int argc, char **argv) {
  const string rs (RE_DOT_STAR RE_NONWHITE_CHAR_CLASS RE_DOT_STAR);
  //  cerr << "Testing regex " << rs << endl;
  const regex nonwhite_re (rs, regex_constants::basic);
  
  smatch sm;
  for (char c = 33; c <= 126; ++c) {
    string s = " " + string(1,c) + " ";
    TEST2 (regex_match (s, sm, nonwhite_re), string(" with s='") + s + "'");
  }

  string w = "   ";
  TEST (!regex_match (w, sm, nonwhite_re));

  const regex stockholm_re (RE_WHITE_OR_EMPTY "#" RE_WHITE_OR_EMPTY "STOCKHOLM" RE_DOT_STAR);
  const regex nexus_re (RE_WHITE_OR_EMPTY "#" RE_WHITE_OR_EMPTY "NEXUS" RE_DOT_STAR);
  const regex fasta_re (RE_WHITE_OR_EMPTY ">" RE_DOT_STAR);
  const regex newick_re (RE_WHITE_OR_EMPTY "\\x28" RE_DOT_STAR);  // '('
  const regex json_re (RE_WHITE_OR_EMPTY "\\x7b" RE_DOT_STAR);    // '{'

  TEST (test_regex (stockholm_re, "# STOCKHOLM 1.0", "## GFF"));
  TEST (test_regex (nexus_re, "# NEXUS", "## GFF"));
  TEST (test_regex (fasta_re, ">seqname", "seqname"));
  TEST (test_regex (newick_re, "(a:1,b:2);", ">abc"));
  TEST (test_regex (json_re, "{a:1,b:2}", "(a:1,b:2)"));

  exit (EXIT_SUCCESS);
}
