#ifndef PAIRHMM_INCLUDED
#define PAIRHMM_INCLUDED

#include "model.h"
#include "logsumexp.h"
#include "profile.h"

struct PairHMM : AlphabetOwner {
  const ProbModel& l;
  const ProbModel& r;
  const LogProbModel logl, logr;
  vector<vector<LogProb> > logRoot;  // log(cptWeight) factored in

  typedef enum { IMM = 0, IMD = 1, IDM = 2, IMI = 3, IIW = 4,
		 TotalStates = 5,
		 SSS = 0, SSI = 3, SIW = 4,
		 EEE = 5  /* EEE should come after all absorbing states, for correct sorting */
  } State;

  // helper methods
  inline int components() const { return logRoot.size(); }

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

  static const char* stateName (State s, bool xAtStart, bool yAtStart);

  // Transition log-probabilities.
  // States {sss,ssi,siw} have same outgoing transition weights as states {imm,imi,iiw}
  // States involving overlapping events (idd,idi) are dropped.
  // Transitions between indistinguishable types of gap (iiw->imd, imi->idm) are also dropped.
  // State iix is dropped because it is only reachable via such an "indistinguishable" transition (imd->iix).
  LogProb imm_imi, imm_iiw, imm_imm, imm_imd, imm_idm, imm_eee;
  LogProb imd_imm, imd_imd, imd_idm, imd_eee;
  LogProb idm_imm, idm_imd, idm_idm, idm_eee;
  LogProb imi_imi, imi_iiw, imi_imm, imi_imd, imi_eee;
  LogProb iiw_iiw, iiw_imm, iiw_idm, iiw_eee;

  // constructor
  PairHMM (const ProbModel& l, const ProbModel& r, const vector<gsl_vector*>& root);

  // helpers
  static vguard<State> states();  // excludes EEE
  static vguard<State> sources (State dest);
  LogProb lpTrans (State src, State dest) const;
};

#endif /* PAIRHMM_INCLUDED */
