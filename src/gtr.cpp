#include "gtr.h"
#include "jsonutil.h"

RateModel gtrModel() {
  RateModel m;
  ParsedJson pj (gtrModelText);
  m.read (pj.value);
  return m;
}

const char* gtrModelText =
"{\n"
"    \"alphabet\": \"ACGT\",\n"
"    \"subrate\" : {\n"
"	\"A\" : { \"C\": 0.3333, \"G\": 0.3333, \"T\": 0.3333 },\n"
"	\"C\" : { \"A\": 0.3333, \"G\": 0.3333, \"T\": 0.3333 },\n"
"	\"G\" : { \"A\": 0.3333, \"C\": 0.3333, \"T\": 0.3333 },\n"
"	\"T\" : { \"A\": 0.3333, \"C\": 0.3333, \"G\": 0.3333 }\n"
"    },\n"
"    \"insrate\" : 0.01,\n"
"    \"delrate\" : 0.01,\n"
"    \"insextprob\" : 0.66,\n"
"    \"delextprob\" : 0.66\n"
"}\n"
;
