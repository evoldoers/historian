#include <cmath>
#include "pairhmm.h"
#include "util.h"

PairHMM::PairHMM (const ProbModel& l, const ProbModel& r, const vector<gsl_vector*>& root)
  : AlphabetOwner (l),
    l (l),
    r (r),
    logl (l),
    logr (r),
    logRoot (log_vector_gsl_vector (root))
{
  for (int cpt = 0; cpt < l.components(); ++cpt)
    for (auto& lr: logRoot[cpt])
      lr += logl.logCptWeight[cpt];

  imm_imi = log (rIns());
  imm_iiw = log (lIns() * rNoIns());
  imm_imm = log (lNoIns() * rNoIns() * lNoDel() * rNoDel());
  imm_imd = log (lNoIns() * rNoIns() * lNoDel() * rDel());
  imm_idm = log (lNoIns() * rNoIns() * lDel() * rNoDel());
  imm_eee = log (lNoIns() * rNoIns());

  imd_imm = log (lNoIns() * lNoDel() * rNoDelExt());
  imd_imd = log (lNoIns() * lNoDel() * rDelExt());
  imd_idm = log (lNoIns() * lDel() * rNoDelExt());
  imd_eee = log (lNoIns() * rNoDelExt());

  idm_imm = log (rNoIns() * lNoDelExt() * rNoDel());
  idm_imd = log (rNoIns() * lNoDelExt() * rDel());
  idm_idm = log (rNoIns() * lDelExt() * rNoDel());
  idm_eee = log (rNoIns() * lNoDelExt());

  imi_imi = log (rInsExt());
  imi_iiw = log (lIns() * rNoInsExt());
  imi_imm = log (lNoIns() * rNoInsExt() * lNoDel() * rNoDel());
  imi_imd = log (lNoIns() * rNoInsExt() * lNoDel() * rDel());
  imi_eee = log (lNoIns() * rNoInsExt());

  iiw_iiw = log (lInsExt());
  iiw_imm = log (lNoInsExt() * lNoDel() * rNoDel());
  iiw_idm = log (lNoInsExt() * lDel() * rNoDel());
  iiw_eee = log (lNoInsExt());
}

LogProb PairHMM::lpTrans (State src, State dest) const {
  switch (src) {
  case IMM:
    switch (dest) {
    case IMM: return imm_imm;
    case IMD: return imm_imd;
    case IDM: return imm_idm;
    case IMI: return imm_imi;
    case IIW: return imm_iiw;
    case EEE: return imm_eee;
    default:
      break;
    }
    break;

  case IMD:
    switch (dest) {
    case IMM: return imd_imm;
    case IMD: return imd_imd;
    case IDM: return imd_idm;
    case EEE: return imd_eee;
    default:
      break;
    }
    break;

  case IDM:
    switch (dest) {
    case IMM: return idm_imm;
    case IMD: return idm_imd;
    case IDM: return idm_idm;
    case EEE: return idm_eee;
    default:
      break;
    }
    break;

  case IMI:
    switch (dest) {
    case IMM: return imi_imm;
    case IMD: return imi_imd;
    case IMI: return imi_imi;
    case IIW: return imi_iiw;
    case EEE: return imi_eee;
    default:
      break;
    }
    break;

  case IIW:
    switch (dest) {
    case IMM: return iiw_imm;
    case IIW: return iiw_iiw;
    case IDM: return iiw_idm;
    case EEE: return iiw_eee;
    default:
      break;
    }
    break;

  default:
    break;
  }
  return -numeric_limits<double>::infinity();
}

vguard<PairHMM::State> PairHMM::states() {
  vguard<State> s = { IMM, IMD, IDM, IMI, IIW };
  return s;
}

vguard<PairHMM::State> PairHMM::sources (State dest) {
  vguard<State> s;
  switch (dest) {
  case IMM:
  case EEE:
    s = { IMM, IMD, IDM, IMI, IIW };
    break;
  case IMD:
    s = { IMM, IMD, IDM, IMI };
    break;
  case IDM:
    s = { IMM, IMD, IDM, IIW };
    break;
  case IMI:
    s = { IMM, IMI };
    break;
  case IIW:
    s = { IMM, IIW, IMI };
    break;
  default:
    break;
  }
  return s;
}

const char* PairHMM::stateName (State s, bool xAtStart, bool yAtStart) {
  switch (s) {
  case IMM: return xAtStart && yAtStart ? "SSS" : "IMM"; break;
  case IMD: return "IMD"; break;
  case IDM: return "IDM"; break;
  case IMI: return xAtStart ? "SSI" : "IMI"; break;
  case IIW: return yAtStart ? "SIW" : "IIW"; break;
  case EEE: return "EEE"; break;
  default: break;
  }
  Abort ("Don't know name of state %u", s);
  return "?";
}

void PairHMM::write (ostream& out) const {
  for (int src = 0; src < (int) TotalStates; ++src)
    for (int dest = 0; dest <= (int) TotalStates; ++dest)
      out << "log P(" << stateName((State)src,false,false) << "->" << stateName((State)dest,false,false) << ") = " << lpTrans((State)src,(State)dest) << endl;
}
