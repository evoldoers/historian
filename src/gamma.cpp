#include <gsl/gsl_cdf.h>
#include "gamma.h"
#include "logger.h"

RateModel makeDiscretizedGammaModel (const RateModel& model, int bins, double shape) {
  Assert (model.components() == 1, "Can't make a discretized gamma model from a model that is already a mixture model");

  RateModel gammaModel (model.alphabet, bins);
  gammaModel.copyIndelParams (model);
  
  vguard<double> rateMultiplier (bins);
  double total = 0;
  for (int c = 0; c < bins; ++c) {
    const double q = ((double) (c + 1)) / ((double) (bins + 1));
    rateMultiplier[c] = gsl_cdf_gamma_Pinv (q, shape, shape);
    total += rateMultiplier[c];
  }
  // adjust so mean rate multiplier = 1
  const double mean = total / (double) bins;
  for (auto& r: rateMultiplier)
    r /= mean;

  LogThisAt(5,"Discretized-gamma rate multipliers: (" << to_string_join(rateMultiplier) << ")" << endl);
  
  for (int c = 0; c < bins; ++c) {
    CheckGsl (gsl_vector_memcpy (gammaModel.insProb[c], model.insProb[0]));
    CheckGsl (gsl_matrix_memcpy (gammaModel.subRate[c], model.subRate[0]));
    CheckGsl (gsl_matrix_scale (gammaModel.subRate[c], rateMultiplier[c]));
  }

  return gammaModel;
}
