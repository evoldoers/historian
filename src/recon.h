#ifndef RECON_INCLUDED
#define RECON_INCLUDED

#include <deque>
#include "knhx.h"
#include "alignpath.h"
#include "model.h"

struct Reconstructor {
  string treeFilename, seqsFilename, modelFilename;
  size_t profileSamples, profileNodeLimit;

  RateModel model;
  ktree_t* tree;
  vguard<FastSeq> seqs;
  
  Reconstructor();
  ~Reconstructor();

  bool parseReconArgs (deque<string>& argvec);

  Alignment reconstruct();
  Alignment loadFilesAndReconstruct();
};


#endif /* PROGALIGN_INCLUDED */
