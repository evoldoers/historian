#include <iostream>
#include <fstream>
#include <string.h>
#include "../src/regexmacros.h"

using namespace std;

int main (int argc, char **argv) {
  const string rs (RE_DOT_STAR RE_NONWHITE_CHAR_CLASS RE_DOT_STAR);
  //  cerr << "Testing regex " << rs << endl;
  const regex nonwhite_re (rs, regex_constants::basic);
  
  smatch sm;
  for (char c = 33; c <= 126; ++c) {
    string s = " " + string(1,c) + " ";
    if (!regex_match (s, sm, nonwhite_re))
      exit (EXIT_FAILURE);
  }

  string w = "   ";
  if (regex_match (w, sm, nonwhite_re))
    exit (EXIT_FAILURE);
  
  exit (EXIT_SUCCESS);
}
