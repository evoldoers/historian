#include <cmath>
#include "pairhmm.h"

PairHMM::PairHMM (const ProbModel& l, const ProbModel& r)
  : l (l),
    r (r)
{
  imm_imi = log (rIns());
  imm_iiw = log (lIns() * rNoIns());
  imm_imm = log (lNoIns() * rNoIns() * lNoDel() * rNoDel());
  imm_imd = log (lNoIns() * rNoIns() * lNoDel() * rDel());
  imm_idm = log (lNoIns() * rNoIns() * lDel() * rNoDel());
  imm_eee = log (lNoIns() * rNoIns());

  imd_iix = log (lIns());
  imd_imm = log (lNoIns() * lNoDel() * rNoDelExt());
  imd_imd = log (lNoIns() * lNoDel() * rDelExt());
  imd_idm = log (lNoIns() * lDel() * rNoDelExt());
  imd_eee = log (lNoIns());

  idm_idi = log (rIns());
  idm_imm = log (rNoIns() * lNoDelExt() * rNoDel());
  idm_imd = log (rNoIns() * lNoDelExt() * rDel());
  idm_idm = log (rNoIns() * lDelExt() * rNoDel());
  idm_eee = log (rNoIns());

  imi_imi = log (rInsExt());
  imi_iiw = log (lIns() * rNoInsExt());
  imi_imm = log (lNoIns() * rNoInsExt() * lNoDel() * rNoDel());
  imi_imd = log (lNoIns() * rNoInsExt() * lNoDel() * rDel());
  imi_idm = log (lNoIns() * rNoInsExt() * lDel() * rNoDel());
  imi_eee = log (lNoIns() * rNoInsExt());

  iiw_iiw = log (lInsExt());
  iiw_imm = log (lNoInsExt() * lNoDel() * rNoDel());
  iiw_imd = log (lNoInsExt() * lNoDel() * rDel());
  iiw_idm = log (lNoInsExt() * lDel() * rNoDel());
  iiw_eee = log (lNoInsExt());

  idi_idi = log (rInsExt());
  idi_imm = log (rNoInsExt() * lNoDelExt() * rNoDel());
  idi_imd = log (rNoInsExt() * lNoDelExt() * rDel());
  idi_idm = log (rNoInsExt() * lDelExt() * rNoDel());
  idi_eee = log (rNoInsExt());

  iix_iix = log (lInsExt());
  iix_imm = log (lNoInsExt() * lNoDel() * rNoDelExt());
  iix_imd = log (lNoInsExt() * lNoDel() * rDelExt());
  iix_idm = log (lNoInsExt() * lDel() * rNoDelExt());
  iix_eee = log (lNoInsExt());
}
