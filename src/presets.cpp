#include "presets.h"
#include "util.h"

#include "ECMrest.h"
#include "ECMunrest.h"
#include "jc.h"
#include "jcrna.h"
#include "lg.h"
#include "wag.h"
#include "jones.h"
#include "dayhoff.h"

RateModel namedModel (const string& n) {
  const string name = tolower(n);
  if (name == "ecmrest")
    return ECMrestModel();
  else if (name == "ecmunrest")
    return ECMunrestModel();
  else if (name == "jc")
    return jcModel();
  else if (name == "jcrna")
    return jcrnaModel();
  else if (name == "lg")
    return lgModel();
  else if (name == "wag")
    return wagModel();
  else if (name == "jtt")
    return jonesModel();
  else if (name == "dayhoff")
    return dayhoffModel();
  else
    Fail ("Unknown model: %s", name.c_str());
  return RateModel();
}
