#include "jcrna.h"
#include "jsonutil.h"

RateModel jcrnaModel() {
  RateModel m;
  ParsedJson pj (jcrnaModelText);
  m.read (pj.value);
  return m;
}

const char* jcrnaModelText =
"{\n"
"    \"alphabet\": \"acgu\",\n"
"    \"wildcard\": \"n\",\n"
"    \"subrate\" : {\n"
"	\"a\" : { \"c\": 0.3333, \"g\": 0.3333, \"u\": 0.3333 },\n"
"	\"c\" : { \"a\": 0.3333, \"g\": 0.3333, \"u\": 0.3333 },\n"
"	\"g\" : { \"a\": 0.3333, \"c\": 0.3333, \"u\": 0.3333 },\n"
"	\"u\" : { \"a\": 0.3333, \"c\": 0.3333, \"g\": 0.3333 }\n"
"    },\n"
"    \"insrate\" : 0.01,\n"
"    \"delrate\" : 0.01,\n"
"    \"insextprob\" : 0.66,\n"
"    \"delextprob\" : 0.66\n"
"}\n"
;
