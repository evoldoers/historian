#ifndef PROGALIGN_INCLUDED
#define PROGALIGN_INCLUDED

#include "knhx.h"
#include "alignpath.h"
#include "model.h"

struct ProgAlignParams {
  size_t profileSamples, profileNodeLimit;

  ProgAlignParams();
  AlignPath reconstructAlignment (ktree_t* tree, const vguard<FastSeq>& seqs, const RateModel& model);
};


#endif /* PROGALIGN_INCLUDED */
