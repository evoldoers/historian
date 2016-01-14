#ifndef PAIRHMM_INCLUDED
#define PAIRHMM_INCLUDED

#include "model.h"

struct PairHMM {
  const RootModel& root;
  const ProbModel& l;
  const ProbModel& r;

  // helper methods
  inline double emit() const { return root.extend; }
  inline double end() const { return 1 - root.extend; }

  inline double lIns() const { return l.ins; }
  inline double lDel() const { return l.del; }
  inline double lInsExt() const { return l.insExt; }
  inline double lDelExt() const { return l.delExt; }

  inline double lNoIns() const { return 1 - l.ins; }
  inline double lNoDel() const { return 1 - l.del; }
  inline double lNoInsExt() const { return 1 - l.insExt; }
  inline double lNoDelExt() const { return 1 - l.delExt; }

  inline double rIns() const { return r.ins; }
  inline double rDel() const { return r.del; }
  inline double rInsExt() const { return r.insExt; }
  inline double rDelExt() const { return r.delExt; }

  inline double rNoIns() const { return 1 - r.ins; }
  inline double rNoDel() const { return 1 - r.del; }
  inline double rNoInsExt() const { return 1 - r.insExt; }
  inline double rNoDelExt() const { return 1 - r.delExt; }

  // transition probabilities
  // States {sss,ssi,siw} have same outgoing transition weights as states {imm,imi,iiw}
  // State imm is dropped.
  double imm_imi, imm_iiw, imm_imm, imm_imd, imm_idm, imm_eee;
  double imd_iix, imd_imm, imd_imd, imd_idm, imd_eee, idm_idi;
  double idm_imm, idm_imd, idm_idm, idm_eee;
  double imi_imi, imi_iiw, imi_imm, imi_imd, imi_idm, imi_eee;
  double iiw_iiw, iiw_imm, iiw_imd, iiw_idm, iiw_eee;
  double idi_idi, idi_imm, idi_imd, idi_idm, idi_eee;
  double iix_iix, iix_imm, iix_imd, iix_idm, iix_eee;

  // constructor
  PairHMM (const RootModel& root, const ProbModel& l, const ProbModel& r);
};

#endif /* PAIRHMM_INCLUDED */
