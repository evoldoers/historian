#ifndef PAIRHMM_INCLUDED
#define PAIRHMM_INCLUDED

#include "model.h"
#include "logsumexp.h"

struct PairHMM : AlphabetOwner {
  const ProbModel& l;
  const ProbModel& r;

  typedef enum { IMM = 0, IMD = 1, IDM = 2, IMI = 3, IDI = 4, IIW = 5, IIX = 6, TotalStates = 7 } State;

  // helper methods
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

  // Transition log-probabilities.
  // States {sss,ssi,siw} have same outgoing transition weights as states {imm,imi,iiw}
  // State idd is dropped.
  LogProb imm_imi, imm_iiw, imm_imm, imm_imd, imm_idm, imm_eee;
  LogProb imd_iix, imd_imm, imd_imd, imd_idm, imd_eee;
  LogProb idm_idi, idm_imm, idm_imd, idm_idm, idm_eee;
  LogProb imi_imi, imi_iiw, imi_imm, imi_imd, imi_idm, imi_eee;
  LogProb iiw_iiw, iiw_imm, iiw_imd, iiw_idm, iiw_eee;
  LogProb idi_idi, idi_imm, idi_imd, idi_idm, idi_eee;
  LogProb iix_iix, iix_imm, iix_imd, iix_idm, iix_eee;

  // constructor
  PairHMM (const ProbModel& l, const ProbModel& r);
};

#endif /* PAIRHMM_INCLUDED */
