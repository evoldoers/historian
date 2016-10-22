#ifndef SIMULATOR_INCLUDED
#define SIMULATOR_INCLUDED

#include <random>

#include "model.h"
#include "tree.h"
#include "alignpath.h"
#include "stockholm.h"
#include "forward.h"

#define SimulatorComponentTag "CPT"

struct Simulator {
  typedef DPMatrix::random_engine random_engine;

  static AlignPath simulateGapsByGillespie (random_engine& generator, const RateModel& model, SeqIdx parentLength, double time, AlignRowIndex parentRowIndex, AlignRowIndex childRowIndex);
  static vguard<FastSeq> simulateSubsByMatExp (random_engine& generator, const RateModel& model, const Tree& tree, const AlignPath& path);  // returns components in qual scores
  static Stockholm simulateTree (random_engine& generator, const RateModel& model, const Tree& tree, SeqIdx rootLength);
};

#endif /* SIMULATOR_INCLUDED */
