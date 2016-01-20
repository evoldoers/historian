#include "forward.h"
#include "util.h"

ForwardMatrix::ForwardMatrix (const Profile& x, const Profile& y, const PairHMM& hmm)
  : x(x),
    y(y),
    hmm(hmm),
    alphSize (hmm.alphabetSize()),
    xSize (x.size()),
    ySize (y.size()),
    subx (x.leftMultiply (hmm.l.subMat)),
    suby (y.leftMultiply (hmm.r.subMat)),
    cellStorage (x.size() * y.size() * PairHMM::TotalStates, -numeric_limits<double>::infinity()),
    absorbScratch (hmm.alphabetSize()),
    insx (x.size(), -numeric_limits<double>::infinity()),
    insy (y.size(), -numeric_limits<double>::infinity())
{
  for (ProfileStateIndex i = 1; i < xSize - 1; ++i)
    insx[i] = logInnerProduct (hmm.l.insVec, x.state[i].lpAbsorb);

  for (ProfileStateIndex j = 1; j < ySize - 1; ++j)
    insy[j] = logInnerProduct (hmm.r.insVec, y.state[j].lpAbsorb);

  cell(0,0,PairHMM::IMM) = 0;
  for (ProfileStateIndex i = 0; i < xSize - 1; ++i) {
    const ProfileState& xState = x.state[i];
    for (ProfileStateIndex j = 0; j < ySize - 1; ++j) {
      const ProfileState& yState = y.state[j];
      if (i > 0) {
	// x-absorbing transitions into IMD, IIW, IIX
	double imd = -numeric_limits<double>::infinity();
	double iiw = -numeric_limits<double>::infinity();
	double iix = -numeric_limits<double>::infinity();
	for (auto xt : xState.in) {
	  const ProfileTransition& xTrans = x.trans[xt];

	  log_accum_exp (imd,
			 log_sum_exp (cell(xTrans.src,j,PairHMM::IMM) + hmm.imm_imd,
				      cell(xTrans.src,j,PairHMM::IMD) + hmm.imd_imd,
				      cell(xTrans.src,j,PairHMM::IDM) + hmm.idm_imd,
				      cell(xTrans.src,j,PairHMM::IMI) + hmm.imi_imd,
				      cell(xTrans.src,j,PairHMM::IIW) + hmm.iiw_imd,
				      cell(xTrans.src,j,PairHMM::IDI) + hmm.idi_imd,
				      cell(xTrans.src,j,PairHMM::IIX) + hmm.iix_imd)
			 + xTrans.lpTrans);

	  log_accum_exp (iiw,
			 log_sum_exp (cell(xTrans.src,j,PairHMM::IMM) + hmm.imm_iiw,
				      cell(xTrans.src,j,PairHMM::IMI) + hmm.imi_iiw,
				      cell(xTrans.src,j,PairHMM::IIW) + hmm.iiw_iiw)
			 + xTrans.lpTrans);

	  log_accum_exp (iix,
			 log_sum_exp (cell(xTrans.src,j,PairHMM::IMD) + hmm.imd_iix,
				      cell(xTrans.src,j,PairHMM::IIX) + hmm.iix_iix)
			 + xTrans.lpTrans);
	}

	cell(i,j,PairHMM::IMD) = imd + insx[i];
	cell(i,j,PairHMM::IIW) = iiw + insx[i];
	cell(i,j,PairHMM::IIX) = iix + insx[i];
      }

      if (j > 0) {
	// y-absorbing transitions into IDM, IMI, IDI
	double idm = -numeric_limits<double>::infinity();
	double imi = -numeric_limits<double>::infinity();
	double idi = -numeric_limits<double>::infinity();
	for (auto yt : yState.in) {
	  const ProfileTransition& yTrans = y.trans[yt];

	  log_accum_exp (idm,
			 log_sum_exp (cell(i,yTrans.src,PairHMM::IMM) + hmm.imm_idm,
				      cell(i,yTrans.src,PairHMM::IMD) + hmm.imd_idm,
				      cell(i,yTrans.src,PairHMM::IDM) + hmm.idm_idm,
				      cell(i,yTrans.src,PairHMM::IMI) + hmm.imi_idm,
				      cell(i,yTrans.src,PairHMM::IIW) + hmm.iiw_idm,
				      cell(i,yTrans.src,PairHMM::IDI) + hmm.idi_idm,
				      cell(i,yTrans.src,PairHMM::IIX) + hmm.iix_idm)
			 + yTrans.lpTrans);

	  log_accum_exp (imi,
			 log_sum_exp (cell(i,yTrans.src,PairHMM::IMM) + hmm.imm_imi,
				      cell(i,yTrans.src,PairHMM::IMI) + hmm.imi_imi)
			 + yTrans.lpTrans);

	  log_accum_exp (idi,
			 log_sum_exp (cell(i,yTrans.src,PairHMM::IDM) + hmm.idm_idi,
				      cell(i,yTrans.src,PairHMM::IDI) + hmm.idi_idi)
			 + yTrans.lpTrans);
	}

	cell(i,j,PairHMM::IDM) = idm + insy[j];
	cell(i,j,PairHMM::IMI) = imi + insy[j];
	cell(i,j,PairHMM::IDI) = idi + insy[j];
      }

      if (i > 0 && j > 0) {
	// xy-absorbing transitions into IMM
	double imm = -numeric_limits<double>::infinity();
	for (auto xt : xState.in) {
	  const ProfileTransition& xTrans = x.trans[xt];
	  for (auto yt : yState.in) {
	    const ProfileTransition& yTrans = y.trans[yt];

	    log_accum_exp (imm,
			   log_sum_exp (cell(xTrans.src,yTrans.src,PairHMM::IMM) + hmm.imm_imm,
					cell(xTrans.src,yTrans.src,PairHMM::IMD) + hmm.imd_imm,
					cell(xTrans.src,yTrans.src,PairHMM::IDM) + hmm.idm_imm,
					cell(xTrans.src,yTrans.src,PairHMM::IMI) + hmm.imi_imm,
					cell(xTrans.src,yTrans.src,PairHMM::IIW) + hmm.iiw_imm,
					cell(xTrans.src,yTrans.src,PairHMM::IDI) + hmm.idi_imm,
					cell(xTrans.src,yTrans.src,PairHMM::IIX) + hmm.iix_imm)
			   + xTrans.lpTrans
			   + yTrans.lpTrans);
	  }
	}

	cell(i,j,PairHMM::IMM) = imm + computeLogProbAbsorb(i,j);
      }
    }
  }

  // transitions into EEE
  for (auto xt : x.end().in) {
    const ProfileTransition& xTrans = x.trans[xt];
    for (auto yt : y.end().in) {
      const ProfileTransition& yTrans = y.trans[yt];
      lpEnd = log_sum_exp (cell(xTrans.src,yTrans.src,PairHMM::IMM) + hmm.imm_eee,
			   cell(xTrans.src,yTrans.src,PairHMM::IMD) + hmm.imd_eee,
			   cell(xTrans.src,yTrans.src,PairHMM::IDM) + hmm.idm_eee,
			   cell(xTrans.src,yTrans.src,PairHMM::IMI) + hmm.imi_eee,
			   cell(xTrans.src,yTrans.src,PairHMM::IIW) + hmm.iiw_eee,
			   cell(xTrans.src,yTrans.src,PairHMM::IDI) + hmm.idi_eee,
			   cell(xTrans.src,yTrans.src,PairHMM::IIX) + hmm.iix_eee)
	+ xTrans.lpTrans
	+ yTrans.lpTrans;
    }
  }
}

ForwardMatrix::CellCoords ForwardMatrix::sampleCell (const map<CellCoords,LogProb>& cellLogProb, random_engine& generator) {
  double ptot = 0;
  for (auto& iter : cellLogProb)
    ptot += exp (iter.second);
  uniform_real_distribution<double> dist (0, ptot);
  double p = dist (generator);
  for (auto& iter : cellLogProb)
    if ((p -= iter.second) <= 0)
      return iter.first;
  Abort ("%s fail",__func__);
  return CellCoords (0, 0, PairHMM::EEE);
}

ForwardMatrix::Path ForwardMatrix::sampleTrace (random_engine& generator) {
  Path path;
  path.push_back (CellCoords (xSize - 1, ySize - 1, PairHMM::EEE));

  const auto states = hmm.states();
  map<CellCoords,LogProb> clp;
  for (auto xt : x.end().in) {
    const ProfileTransition& xTrans = x.trans[xt];
    for (auto yt : y.end().in) {
      const ProfileTransition& yTrans = y.trans[yt];
      for (auto s : states)
	clp[CellCoords(xTrans.src,yTrans.src,s)] = cell(xTrans.src,yTrans.src,s) + hmm.lpTrans (s, PairHMM::EEE) + xTrans.lpTrans + yTrans.lpTrans;
    }
  }

  CellCoords current;
  while (true) {
    current = sampleCell (clp, generator);
    path.push_front (current);
    if (current.xpos == 0 && current.ypos == 0)
      break;

    clp.clear();
    switch (current.state) {
      // x-absorbing transitions into IMD, IIW, IIX
    case PairHMM::IMD:
    case PairHMM::IIW:
    case PairHMM::IIX:
      for (auto xt : x.state[current.xpos].in)
	for (auto s : states)
	  clp[CellCoords(x.trans[xt].src,current.ypos,s)] = cell(x.trans[xt].src,current.ypos,s) + hmm.lpTrans(s,current.state) + x.trans[xt].lpTrans;
      break;

      // y-absorbing transitions into IDM, IMI, IDI
    case PairHMM::IDM:
    case PairHMM::IMI:
    case PairHMM::IDI:
      for (auto yt : y.state[current.ypos].in)
	for (auto s : states)
	  clp[CellCoords(current.xpos,y.trans[yt].src,s)] = cell(current.xpos,y.trans[yt].src,s) + hmm.lpTrans(s,current.state) + y.trans[yt].lpTrans;
      break;

      // xy-absorbing transitions into IMM
    case PairHMM::IMM:
      for (auto xt : x.state[current.xpos].in)
	for (auto yt : y.state[current.ypos].in)
	  for (auto s : states)
	    clp[CellCoords(x.trans[xt].src,y.trans[yt].src,s)] = cell(x.trans[xt].src,y.trans[yt].src,s) + hmm.lpTrans(s,current.state) + x.trans[xt].lpTrans + y.trans[yt].lpTrans;
      break;

    default:
      Abort ("%s fail",__func__);
      break;
    }
  }
  
  return path;
}
