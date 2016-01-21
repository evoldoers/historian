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
    insy (y.size(), -numeric_limits<double>::infinity()),
    delx (x.size(), -numeric_limits<double>::infinity()),
    dely (y.size(), -numeric_limits<double>::infinity())
{
  for (ProfileStateIndex i = 1; i < xSize - 1; ++i) {
    insx[i] = logInnerProduct (hmm.l.insVec, x.state[i].lpAbsorb);
    delx[i] = logInnerProduct (hmm.root, subx.state[i].lpAbsorb);
  }

  for (ProfileStateIndex j = 1; j < ySize - 1; ++j) {
    insy[j] = logInnerProduct (hmm.r.insVec, y.state[j].lpAbsorb);
    dely[j] = logInnerProduct (hmm.root, suby.state[j].lpAbsorb);
  }

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

	cell(i,j,PairHMM::IMD) = imd + delx[i];
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

	cell(i,j,PairHMM::IDM) = idm + dely[j];
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
    clp = sourceCells (current);
  }
  
  return path;
}

map<ForwardMatrix::CellCoords,LogProb> ForwardMatrix::sourceCells (const CellCoords& destCell) {
  map<CellCoords,LogProb> clp;
  switch (destCell.state) {
    // x-absorbing transitions into IMD, IIW, IIX
  case PairHMM::IMD:
  case PairHMM::IIW:
  case PairHMM::IIX:
    for (auto xt : x.state[destCell.xpos].in)
      for (auto s : hmm.states())
	clp[CellCoords(x.trans[xt].src,destCell.ypos,s)] = cell(x.trans[xt].src,destCell.ypos,s) + hmm.lpTrans(s,destCell.state) + x.trans[xt].lpTrans;
    break;

    // y-absorbing transitions into IDM, IMI, IDI
  case PairHMM::IDM:
  case PairHMM::IMI:
  case PairHMM::IDI:
    for (auto yt : y.state[destCell.ypos].in)
      for (auto s : hmm.states())
	clp[CellCoords(destCell.xpos,y.trans[yt].src,s)] = cell(destCell.xpos,y.trans[yt].src,s) + hmm.lpTrans(s,destCell.state) + y.trans[yt].lpTrans;
    break;

    // xy-absorbing transitions into IMM
  case PairHMM::IMM:
    for (auto xt : x.state[destCell.xpos].in)
      for (auto yt : y.state[destCell.ypos].in)
	for (auto s : hmm.states())
	  clp[CellCoords(x.trans[xt].src,y.trans[yt].src,s)] = cell(x.trans[xt].src,y.trans[yt].src,s) + hmm.lpTrans(s,destCell.state) + x.trans[xt].lpTrans + y.trans[yt].lpTrans;
    break;

  default:
    Abort ("%s fail",__func__);
    break;
  }
  
  return clp;
}

bool ForwardMatrix::CellCoords::isAbsorbing() const {
  return state == PairHMM::IMM || state == PairHMM::IMD || state == PairHMM::IDM;
}

LogProb ForwardMatrix::eliminatedLogProbAbsorb (const CellCoords& cell) const {
  switch (cell.state) {
  case PairHMM::IIW:
  case PairHMM::IIX:
    return insx[cell.xpos];
    break;

  case PairHMM::IMI:
  case PairHMM::IDI:
    return insy[cell.ypos];
    break;

  case PairHMM::IMM:
  case PairHMM::IMD:
  case PairHMM::IDM:
  default:
    Abort ("%s fail",__func__);
    break;
  }

  return 0;
}

Profile ForwardMatrix::makeProfile (const set<CellCoords>& cells, AlignRowIndex rowIndex) {
  Profile prof;

  const CellCoords startCell (0, 0, PairHMM::SSS), endCell (CellCoords (xSize - 1, ySize - 1, PairHMM::EEE));
  Assert (cells.find (startCell) != cells.end(), "Missing SSS");
  Assert (cells.find (endCell) != cells.end(), "Missing EEE");

  // build states
  // retain only start, end, and absorbing states
  map<CellCoords,ProfileStateIndex> profStateIndex;
  map<CellCoords,AlignPath> elimAlignPath;
  profStateIndex[startCell] = 0;
  prof.state.push_back (ProfileState());

  for (const auto& c : cells)
    if (c.isAbsorbing()) {
      // cell is to be retained
      profStateIndex[c] = prof.state.size();
      prof.state.push_back (ProfileState());
      switch (c.state) {
      case PairHMM::IMM:
	initAbsorbScratch (c.xpos, c.ypos);
	prof.state.back().lpAbsorb = absorbScratch;
	prof.state.back().alignPath = mergeAlignments (x.state[c.xpos].alignPath, y.state[c.ypos].alignPath);
	break;
      case PairHMM::IMD:
	prof.state.back().lpAbsorb = subx.state[c.xpos].lpAbsorb;
	prof.state.back().alignPath = x.state[c.xpos].alignPath;
	break;
      case PairHMM::IDM:
	prof.state.back().lpAbsorb = suby.state[c.ypos].lpAbsorb;
	prof.state.back().alignPath = y.state[c.ypos].alignPath;
	break;
      default:
	Abort ("%s fail",__func__);
	break;
      }
      prof.state.back().alignPath[rowIndex].push_back (true);
    } else {
      // cell is to be eliminated
      switch (c.state) {
      case PairHMM::IIW:
      case PairHMM::IIX:
	elimAlignPath[c] = x.state[c.xpos].alignPath;
	break;
      case PairHMM::IMI:
      case PairHMM::IDI:
	elimAlignPath[c] = y.state[c.ypos].alignPath;
	break;
      default:
	break;
      }
    }

  profStateIndex[endCell] = prof.state.size();
  prof.state.push_back (ProfileState());

  // build transitions
  map<CellCoords,map<ProfileStateIndex,ProfileTransition> > effTrans;  // includes both direct transitions to retained states, and transitions representing sum-over-eliminated state paths to retained states
  for (auto iter = cells.crbegin(); iter != cells.crend(); ++iter) {
    const CellCoords& cell = *iter;
    const map<CellCoords,LogProb>& slp = sourceCells (cell);
    if (profStateIndex.find(cell) != profStateIndex.end()) {
      // cell is to be retained
      const ProfileStateIndex idx = profStateIndex[cell];
      vector<ProfileTransitionIndex>& out = prof.state[profStateIndex[cell]].out;
      for (const auto& effTransIter : effTrans[cell]) {
	out.push_back (prof.trans.size());
	prof.trans.push_back (effTransIter.second);
	prof.trans.back().src = idx;
      }
      for (auto slpIter : slp) {
	ProfileTransition trans;
	trans.dest = idx;
	trans.lpTrans = slpIter.second;
	effTrans[slpIter.first][idx] = trans;
      }
    } else {
      // cell is to be eliminated
      const auto& cellEffTrans = effTrans[cell];
      const AlignPath& cellAlignPath = elimAlignPath[cell];
      const LogProb cellLogProbAbsorb = eliminatedLogProbAbsorb (cell);
      for (auto slpIter : slp) {
	const CellCoords& srcCell = slpIter.first;
	const LogProb srcCellLogProbTrans = slpIter.second;
	auto& srcEffTrans = effTrans[srcCell];
	for (auto cellEffTransIter : cellEffTrans) {
	  const ProfileStateIndex& destIdx = cellEffTransIter.first;
	  const ProfileTransition& cellDestTrans = cellEffTransIter.second;
	  if (srcEffTrans.find(destIdx) == srcEffTrans.end()) {
	    // TODO: track best alignPath, instead of just using first one encountered
	    ProfileTransition trans;
	    trans.dest = destIdx;
	    trans.lpTrans = -numeric_limits<double>::infinity();
	    trans.alignPath = concatenateAlignments (cellAlignPath, cellDestTrans.alignPath);
	    srcEffTrans[destIdx] = trans;
	  }
	  log_accum_exp (srcEffTrans[destIdx].lpTrans, srcCellLogProbTrans + cellLogProbAbsorb + cellDestTrans.lpTrans);
	}
      }
    }
  }
  
  // WRITE ME
  return prof;
}

