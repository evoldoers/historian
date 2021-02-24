#include "jc.h"
#include "jsonutil.h"

RateModel jcModel() {
  RateModel m;
  ParsedJson pj (jcModelText);
  m.read (pj.value);
  return m;
}

const char* jcModelText =
"{\n"
"    \"alphabet\": \"acgt\",\n"
"    \"wildcard\": \"n\",\n"
"    \"subrate\" : {\n"
"	\"a\" : { \"c\": 0.3333, \"g\": 0.3333, \"t\": 0.3333 },\n"
"	\"c\" : { \"a\": 0.3333, \"g\": 0.3333, \"t\": 0.3333 },\n"
"	\"g\" : { \"a\": 0.3333, \"c\": 0.3333, \"t\": 0.3333 },\n"
"	\"t\" : { \"a\": 0.3333, \"c\": 0.3333, \"g\": 0.3333 }\n"
"    },\n"
"    \"insrate\" : 0.01,\n"
"    \"delrate\" : 0.01,\n"
"    \"insextprob\" : 0.66,\n"
"    \"delextprob\" : 0.66\n"
"}\n"
;
