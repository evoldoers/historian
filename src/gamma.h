#ifndef GAMMA_INCLUDED
#define GAMMA_INCLUDED

#include "model.h"

RateModel makeDiscretizedGammaModel (const RateModel& model, int bins, double shape);

#endif /* GAMMA_INCLUDED */
