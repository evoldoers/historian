#include "model.h"

RateModelBasis::RateModelBasis (const RateModel& model)
  : model (model)
{
  // WRITE ME
}

RootModel RateModelBasis::defaultRootModel() const {
  RootModel r;
  r.init = eqm;
  r.extend = r.end = 1;
  return r;
}

gsl_matrix RateModelBasis::exp (double t) const {
  gsl_matrix m = evec;
  // WRITE ME
  return m;
}

