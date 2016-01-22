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
    dely (y.size(), -numeric_limits<double>::infinity()),
    startCell (0, 0, PairHMM::SSS),
    endCell (xSize - 1, ySize - 1, PairHMM::EEE)
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
	// x-absorbing transitions into IMD, IIW
	double imd = -numeric_limits<double>::infinity();
	double iiw = -numeric_limits<double>::infinity();
	for (auto xt : xState.in) {
	  const ProfileTransition& xTrans = x.trans[xt];

	  log_accum_exp (imd,
			 log_sum_exp (cell(xTrans.src,j,PairHMM::IMM) + hmm.imm_imd,
				      cell(xTrans.src,j,PairHMM::IMD) + hmm.imd_imd,
				      cell(xTrans.src,j,PairHMM::IDM) + hmm.idm_imd,
				      cell(xTrans.src,j,PairHMM::IMI) + hmm.imi_imd)
			 + xTrans.lpTrans);

	  log_accum_exp (iiw,
			 log_sum_exp (cell(xTrans.src,j,PairHMM::IMM) + hmm.imm_iiw,
				      cell(xTrans.src,j,PairHMM::IMI) + hmm.imi_iiw,
				      cell(xTrans.src,j,PairHMM::IIW) + hmm.iiw_iiw)
			 + xTrans.lpTrans);
	}

	cell(i,j,PairHMM::IMD) = imd + delx[i];
	cell(i,j,PairHMM::IIW) = iiw + insx[i];
      }

      if (j > 0) {
	// y-absorbing transitions into IDM, IMI
	double idm = -numeric_limits<double>::infinity();
	double imi = -numeric_limits<double>::infinity();
	for (auto yt : yState.in) {
	  const ProfileTransition& yTrans = y.trans[yt];

	  log_accum_exp (idm,
			 log_sum_exp (cell(i,yTrans.src,PairHMM::IMM) + hmm.imm_idm,
				      cell(i,yTrans.src,PairHMM::IMD) + hmm.imd_idm,
				      cell(i,yTrans.src,PairHMM::IDM) + hmm.idm_idm,
				      cell(i,yTrans.src,PairHMM::IIW) + hmm.iiw_idm)
			 + yTrans.lpTrans);

	  log_accum_exp (imi,
			 log_sum_exp (cell(i,yTrans.src,PairHMM::IMM) + hmm.imm_imi,
				      cell(i,yTrans.src,PairHMM::IMI) + hmm.imi_imi)
			 + yTrans.lpTrans);
	}

	cell(i,j,PairHMM::IDM) = idm + dely[j];
	cell(i,j,PairHMM::IMI) = imi + insy[j];
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
					cell(xTrans.src,yTrans.src,PairHMM::IIW) + hmm.iiw_imm)
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
			   cell(xTrans.src,yTrans.src,PairHMM::IIW) + hmm.iiw_eee)
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
  return CellCoords();
}

ForwardMatrix::CellCoords ForwardMatrix::bestCell (const map<CellCoords,LogProb>& cellLogProb) {
  CellCoords best;
  double pBest = -numeric_limits<double>::infinity();
  for (auto& iter : cellLogProb)
    if (iter.second > pBest) {
      pBest = iter.second;
      best = iter.first;
    }
  return best;
}

ForwardMatrix::Path ForwardMatrix::sampleTrace (random_engine& generator) {
  Path path;
  path.push_back (endCell);

  map<CellCoords,LogProb> clp = sourceCells (endCell);
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

ForwardMatrix::Path ForwardMatrix::bestTrace() {
  Path path;
  path.push_back (endCell);

  map<CellCoords,LogProb> clp = sourceCells (endCell);
  CellCoords current;
  while (true) {
    current = bestCell (clp);
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
    // x-absorbing transitions into IMD, IIW
  case PairHMM::IMD:
  case PairHMM::IIW:
    for (auto xt : x.state[destCell.xpos].in)
      for (auto s : hmm.states())
	clp[CellCoords(x.trans[xt].src,destCell.ypos,s)] = cell(x.trans[xt].src,destCell.ypos,s) + hmm.lpTrans(s,destCell.state) + x.trans[xt].lpTrans;
    break;

    // y-absorbing transitions into IDM, IMI
  case PairHMM::IDM:
  case PairHMM::IMI:
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

    // null transitions into EEE
  case PairHMM::EEE:
    if (destCell.xpos == xSize - 1 && destCell.ypos == ySize - 1)
      for (auto xt : x.end().in)
	for (auto yt : y.end().in)
	  for (auto s : hmm.states())
	    clp[CellCoords(x.trans[xt].src,y.trans[yt].src,s)] = cell(x.trans[xt].src,y.trans[yt].src,s) + hmm.lpTrans (s,destCell.state) + x.trans[xt].lpTrans + y.trans[yt].lpTrans;
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
    return insx[cell.xpos];
    break;

  case PairHMM::IMI:
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

ForwardMatrix::EffectiveTransition::EffectiveTransition()
  : lpPath (-numeric_limits<double>::infinity()),
    lpBestAlignPath (-numeric_limits<double>::infinity())
{ }

AlignPath ForwardMatrix::cellAlignPath (const CellCoords& c, AlignRowIndex rowIndex) const {
  AlignPath alignPath;
  if (c.isAbsorbing()) {
    switch (c.state) {
      alignPath = alignPathUnion (x.state[c.xpos].alignPath, y.state[c.ypos].alignPath);
      break;
    case PairHMM::IMD:
      alignPath = x.state[c.xpos].alignPath;
      break;
    case PairHMM::IDM:
      alignPath = y.state[c.ypos].alignPath;
      break;
    default:
      Abort ("%s fail",__func__);
      break;
    }
    alignPath[rowIndex].push_back (true);
  }
  return alignPath;
}

AlignPath ForwardMatrix::transitionAlignPath (const CellCoords& src, const CellCoords& dest) const {
  AlignPath path;
  if (src.xpos != dest.xpos)
    path = x.getTrans(src.xpos,dest.xpos)->alignPath;
  if (src.ypos != dest.ypos)
    path = alignPathConcat (path, y.getTrans(src.ypos,dest.ypos)->alignPath);
  return path;
}

Profile ForwardMatrix::makeProfile (const set<CellCoords>& cells, AlignRowIndex rowIndex) {
  Profile prof;

  Assert (cells.find (startCell) != cells.end(), "Missing SSS");
  Assert (cells.find (endCell) != cells.end(), "Missing EEE");

  // build states
  // retain only start, end, and absorbing cells
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
	break;
      case PairHMM::IMD:
	prof.state.back().lpAbsorb = subx.state[c.xpos].lpAbsorb;
	break;
      case PairHMM::IDM:
	prof.state.back().lpAbsorb = suby.state[c.ypos].lpAbsorb;
	break;
      default:
	Abort ("%s fail",__func__);
	break;
      }
      prof.state.back().alignPath = cellAlignPath (c, rowIndex);
    } else {
      // cell is to be eliminated
      switch (c.state) {
      case PairHMM::IIW:
	elimAlignPath[c] = x.state[c.xpos].alignPath;
	break;
      case PairHMM::IMI:
	elimAlignPath[c] = y.state[c.ypos].alignPath;
	break;
      default:
	break;
      }
    }

  profStateIndex[endCell] = prof.state.size();
  prof.state.push_back (ProfileState());

  // Calculate log-probabilities of "effective transitions" from cells to retained-cells.
  // Each effective transition represents a sum over paths.
  // A path is either a single direct transition (from source cell to destination retained-cell),
  // or a series of transitions starting from the source cell,
  // passing through one or more eliminated-cells, and stopping at the destination retained-cell.
  map<CellCoords,map<ProfileStateIndex,EffectiveTransition> > effTrans;  // effTrans[srcCell][destStateIdx]
  for (auto iter = cells.crbegin(); iter != cells.crend(); ++iter) {
    const CellCoords& cell = *iter;
    const map<CellCoords,LogProb>& slp = sourceCells (cell);
    if (profStateIndex.find(cell) != profStateIndex.end()) {
      // cell is to be retained. Incoming & outgoing paths can be kept separate
      const ProfileStateIndex cellIdx = profStateIndex[cell];
      for (const auto slpIter : slp) {
	const CellCoords& src = slpIter.first;
	const LogProb srcCellLogProbTrans = slpIter.second;
	EffectiveTransition& eff = effTrans[src][cellIdx];
	eff.lpPath = eff.lpBestAlignPath = srcCellLogProbTrans;
      }
    } else {
      // cell is to be eliminated. Connect incoming transitions & outgoing paths, summing cell out
      const auto& cellEffTrans = effTrans[cell];
      const AlignPath& cellAlignPath = elimAlignPath[cell];
      const LogProb cellLogProbAbsorb = eliminatedLogProbAbsorb (cell);
      for (auto slpIter : slp) {
	const CellCoords& src = slpIter.first;
	const LogProb srcCellLogProbTrans = slpIter.second;
	auto& srcEffTrans = effTrans[src];
	for (const auto cellEffTransIter : cellEffTrans) {
	  const ProfileStateIndex& destIdx = cellEffTransIter.first;
	  const EffectiveTransition& cellDestEffTrans = cellEffTransIter.second;
	  EffectiveTransition& srcDestEffTrans = srcEffTrans[destIdx];
	  log_accum_exp (srcDestEffTrans.lpPath, srcCellLogProbTrans + cellLogProbAbsorb + cellDestEffTrans.lpPath);
	  // we also want to keep the best single alignment consistent w/this set of paths
	  const LogProb srcDestLogProbBestAlignPath = srcCellLogProbTrans + cellLogProbAbsorb + cellDestEffTrans.lpBestAlignPath;
	  if (srcDestLogProbBestAlignPath > srcDestEffTrans.lpBestAlignPath) {
	    srcDestEffTrans.lpBestAlignPath = srcDestLogProbBestAlignPath;
	    srcDestEffTrans.bestAlignPath = alignPathConcat (transitionAlignPath(src,cell), cellAlignPath, cellDestEffTrans.bestAlignPath);
	  }
	}
      }
    }
  }

  // populate outgoing & incoming transitions for each state
  for (const auto profStateIter : profStateIndex) {
    const CellCoords& cell = profStateIter.first;
    const ProfileStateIndex srcIdx = profStateIter.second;
    vector<ProfileTransitionIndex>& srcOut = prof.state[srcIdx].out;
    for (const auto effTransIter : effTrans[cell]) {
      const ProfileStateIndex destIdx = effTransIter.first;
      const EffectiveTransition& srcDestEffTrans = effTransIter.second;
      vector<ProfileTransitionIndex>& destIn = prof.state[destIdx].in;
      const ProfileTransitionIndex transIdx = prof.trans.size();
      ProfileTransition trans;
      trans.src = srcIdx;
      trans.dest = destIdx;
      trans.lpTrans = srcDestEffTrans.lpPath;
      trans.alignPath = srcDestEffTrans.bestAlignPath;
      prof.trans.push_back (trans);
      srcOut.push_back (transIdx);
      destIn.push_back (transIdx);
    }
  }
  
  return prof;
}

Profile ForwardMatrix::sampleProfile (size_t profileSamples, size_t maxCells, random_engine& generator, AlignRowIndex rowIndex) {
  map<CellCoords,size_t> cellCount;
  const Path best = bestTrace();
  for (auto& c : best)
    cellCount[c] = 2;  // avoid dropping these cells
  for (size_t n = 0; n < profileSamples && cellCount.size() < maxCells; ++n) {
    const Path sampled = sampleTrace (generator);
    for (auto& c : sampled)
      ++cellCount[c];
  }
  set<CellCoords> profCells;
  for (const auto& cc : cellCount)
    if (cc.second > 1)
      profCells.insert (cc.first);
  return makeProfile (profCells, rowIndex);
}

Profile ForwardMatrix::bestProfile (AlignRowIndex rowIndex) {
  const Path best = bestTrace();
  const set<CellCoords> profCells (best.begin(), best.end());
  return makeProfile (profCells, rowIndex);
}

