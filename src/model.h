#ifndef MODEL_INCLUDED
#define MODEL_INCLUDED

#include <gsl/gsl_matrix.h>
#include <string>
#include "gason.h"

using namespace std;

struct AlphabetOwner {
  string alphabet;
  size_t alphabetSize() const { return alphabet.size(); }
};

struct RootModel : AlphabetOwner {
  gsl_vector init;
  double extend, end;  // Strict (geometric) normalization: end = 1 - extend. Weak (flat): end = extend = 1
};

struct RateModel : AlphabetOwner {
  double insRate, delRate, insExtProb, delExtProb;
  gsl_matrix subRate;

  void read (const JsonValue& json);
};

struct RateModelBasis {
  const RateModel& model;
  gsl_vector eqm;
  gsl_vector eval;
  gsl_matrix evec;
  RateModelBasis (const RateModel& model);
  RootModel defaultRootModel() const;
  gsl_matrix exp (double t) const;
};

struct ProbModel : AlphabetOwner {
  double ins, del, insExt, delExt;
  gsl_matrix sub;
  ProbModel (const RateModelBasis& basis, double t);
};

#endif /* MODEL_INCLUDED */
