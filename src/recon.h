#ifndef RECON_INCLUDED
#define RECON_INCLUDED

#include "knhx.h"
#include "alignpath.h"
#include "model.h"

struct Reconstructor {
  size_t profileSamples, profileNodeLimit;

  Reconstructor();
  AlignPath reconstruct (ktree_t* tree, const vguard<FastSeq>& seqs, const RateModel& model);
};


#endif /* PROGALIGN_INCLUDED */
