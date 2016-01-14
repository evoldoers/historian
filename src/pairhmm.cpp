#include <cmath>
#include "pairhmm.h"

PairHMM::PairHMM (const RootModel& root, const ProbModel& l, const ProbModel& r)
  : root (root),
    l (l),
    r (r)
{
  imm_imi = log (rIns());
  imm_iiw = log (lIns() * rNoIns());
  imm_imm = log (lNoIns() * rNoIns() * emit() * lNoDel() * rNoDel());
  imm_imd = log (lNoIns() * rNoIns() * emit() * lNoDel() * rDel());
  imm_idm = log (lNoIns() * rNoIns() * emit() * lDel() * rNoDel());
  imm_eee = log (lNoIns() * rNoIns() * end());

  imd_iix = log (lIns());
  imd_imm = log (lNoIns() * emit() * lNoDel() * rNoDelExt());
  imd_imd = log (lNoIns() * emit() * lNoDel() * rDelExt());
  imd_idm = log (lNoIns() * emit() * lDel() * rNoDelExt());
  imd_eee = log (lNoIns() * end());

  idm_idi = log (rIns());
  idm_imm = log (rNoIns() * emit() * lNoDelExt() * rNoDel());
  idm_imd = log (rNoIns() * emit() * lNoDelExt() * rDel());
  idm_idm = log (rNoIns() * emit() * lDelExt() * rNoDel());
  idm_eee = log (rNoIns() * end());

  imi_imi = log (rInsExt());
  imi_iiw = log (lIns() * rNoInsExt());
  imi_imm = log (lNoIns() * rNoInsExt() * emit() * lNoDel() * rNoDel());
  imi_imd = log (lNoIns() * rNoInsExt() * emit() * lNoDel() * rDel());
  imi_idm = log (lNoIns() * rNoInsExt() * emit() * lDel() * rNoDel());
  imi_eee = log (lNoIns() * rNoInsExt() * end());

  iiw_iiw = log (lInsExt());
  iiw_imm = log (lNoInsExt() * emit() * lNoDel() * rNoDel());
  iiw_imd = log (lNoInsExt() * emit() * lNoDel() * rDel());
  iiw_idm = log (lNoInsExt() * emit() * lDel() * rNoDel());
  iiw_eee = log (lNoInsExt() * end());

  idi_idi = log (rInsExt());
  idi_imm = log (rNoInsExt() * emit() * lNoDelExt() * rNoDel());
  idi_imd = log (rNoInsExt() * emit() * lNoDelExt() * rDel());
  idi_idm = log (rNoInsExt() * emit() * lDelExt() * rNoDel());
  idi_eee = log (rNoInsExt() * end());

  iix_iix = log (lInsExt());
  iix_imm = log (lNoInsExt() * emit() * lNoDel() * rNoDelExt());
  iix_imd = log (lNoInsExt() * emit() * lNoDel() * rDelExt());
  iix_idm = log (lNoInsExt() * emit() * lDel() * rNoDelExt());
  iix_eee = log (lNoInsExt() * end());
}
