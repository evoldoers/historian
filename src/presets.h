#ifndef PRESETS_INCLUDED
#define PRESETS_INCLUDED

#include "model.h"

RateModel defaultAminoModel();
RateModel defaultCodonModel();
RateModel defaultBaseModel();

RateModel namedModel (const string& name);

#endif /* PRESETS_INCLUDED */
