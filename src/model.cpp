#include "model.h"

RateModelBasis::RateModelBasis (const RateModel& model)
  : model (model)
{
  // WRITE ME
}

gsl_matrix RateModelBasis::exp (double t) const {
  gsl_matrix m = evec;
  // WRITE ME
  return m;
}

