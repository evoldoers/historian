#include <gsl/gsl_math.h>

#include "forward.h"
#include "util.h"
#include "logger.h"
#include "tree.h"
#include "sumprod.h"

#define FWD_BACK_ERROR_TOLERANCE .01

DPMatrix::DPMatrix (const Profile& x, const Profile& y, const PairHMM& hmm, const GuideAlignmentEnvelope& env)
  : x(x),
    y(y),
    hmm(hmm),
    alphSize ((AlphTok) hmm.alphabetSize()),
    xSize (x.size()),
    ySize (y.size()),
    subx (x.leftMultiply (hmm.l.subMat)),
    suby (y.leftMultiply (hmm.r.subMat)),
    cellStorage (x.size()),
    absorbScratch (hmm.alphabetSize()),
    insx (x.size(), -numeric_limits<double>::infinity()),
    insy (y.size(), -numeric_limits<double>::infinity()),
    rootsubx (x.size(), -numeric_limits<double>::infinity()),
    rootsuby (y.size(), -numeric_limits<double>::infinity()),
    startCell (0, 0, PairHMM::SSS),
    endCell (xSize - 1, ySize - 1, PairHMM::EEE),
    envelope (env),
    xClosestLeafPos (xSize, 0),
    yClosestLeafPos (ySize, 0)
{
  if (env.initialized()) {
    for (ProfileStateIndex i = 1; i < xSize; ++i)
      xClosestLeafPos[i] = x.state[i].seqCoords.at(env.row1);

    for (ProfileStateIndex j = 1; j < ySize; ++j)
      yClosestLeafPos[j] = y.state[j].seqCoords.at(env.row2);
  }

  for (ProfileStateIndex i = 1; i < xSize - 1; ++i)
    if (!x.state[i].isNull()) {
      insx[i] = logInnerProduct (hmm.logl.logInsProb, x.state[i].lpAbsorb);
      rootsubx[i] = logInnerProduct (hmm.logRoot, subx.state[i].lpAbsorb);
    }

  for (ProfileStateIndex j = 1; j < ySize - 1; ++j)
    if (!y.state[j].isNull()) {
      insy[j] = logInnerProduct (hmm.logr.logInsProb, y.state[j].lpAbsorb);
      rootsuby[j] = logInnerProduct (hmm.logRoot, suby.state[j].lpAbsorb);
    }
}

ForwardMatrix::ForwardMatrix (const Profile& x, const Profile& y, const PairHMM& hmm, AlignRowIndex parentRowIndex, const GuideAlignmentEnvelope& env, SumProduct* sumProd)
  : DPMatrix (x, y, hmm, env),
    parentRowIndex (parentRowIndex),
    sumProd (sumProd)
{
  lpStart() = 0;

  ProgressLog (plog, 4);
  plog.initProgress ("Forward algorithm (%s vs %s)", x.name.c_str(), y.name.c_str());

  for (ProfileStateIndex i = 0; i < xSize - 1; ++i) {
    const ProfileState& xState = x.state[i];

    plog.logProgress (i / (double) (xSize - 2), "state %d/%d", i + 1, xSize);

    for (ProfileStateIndex j = 0; j < ySize - 1; ++j) {
      const ProfileState& yState = y.state[j];

      if (envelope.inRange (xClosestLeafPos[i], yClosestLeafPos[j])) {

	XYCell& dest = xyCell(i,j);
	double& imm = dest(PairHMM::IMM);
	double& imd = dest(PairHMM::IMD);
	double& idm = dest(PairHMM::IDM);
	double& imi = dest(PairHMM::IMI);
	double& iiw = dest(PairHMM::IIW);

	if (!xState.isNull()) {
	  // x-absorbing transitions into IMD, IIW
	  for (auto xt : xState.in) {
	    const ProfileTransition& xTrans = x.trans[xt];
	    const XYCell& src = xyCell(xTrans.src,j);

	    log_accum_exp (imd,
			   log_sum_exp (src(PairHMM::IMM) + hmm.imm_imd,
					src(PairHMM::IMD) + hmm.imd_imd,
					src(PairHMM::IDM) + hmm.idm_imd,
					src(PairHMM::IMI) + hmm.imi_imd)
			   + xTrans.lpTrans);

	    log_accum_exp (iiw,
			   log_sum_exp (src(PairHMM::IMM) + hmm.imm_iiw,
					src(PairHMM::IMI) + hmm.imi_iiw,
					src(PairHMM::IIW) + hmm.iiw_iiw)
			   + xTrans.lpTrans);
	  }

	  imd += rootsubx[i];
	  iiw += insx[i];

	} else {
	  // x-nonabsorbing transitions in IMD, IIW
	  for (auto xt : xState.in) {
	    const ProfileTransition& xTrans = x.trans[xt];
	    const XYCell& src = xyCell(xTrans.src,j);
	    log_accum_exp (imd, src(PairHMM::IMD) + xTrans.lpTrans);
	    log_accum_exp (iiw, src(PairHMM::IIW) + xTrans.lpTrans);
	  }
	}

	if (!yState.isNull()) {
	  // y-absorbing transitions into IDM, IMI
	  for (auto yt : yState.in) {
	    const ProfileTransition& yTrans = y.trans[yt];
	    const XYCell& src = xyCell(i,yTrans.src);

	    log_accum_exp (idm,
			   log_sum_exp (src(PairHMM::IMM) + hmm.imm_idm,
					src(PairHMM::IMD) + hmm.imd_idm,
					src(PairHMM::IDM) + hmm.idm_idm,
					src(PairHMM::IIW) + hmm.iiw_idm)
			   + yTrans.lpTrans);

	    log_accum_exp (imi,
			   log_sum_exp (src(PairHMM::IMM) + hmm.imm_imi,
					src(PairHMM::IMI) + hmm.imi_imi)
			   + yTrans.lpTrans);
	  }

	  idm += rootsuby[j];
	  imi += insy[j];

	} else {
	  // y-nonabsorbing transitions in IDM, IMI
	  for (auto yt : yState.in) {
	    const ProfileTransition& yTrans = y.trans[yt];
	    const XYCell& src = xyCell(i,yTrans.src);
	    log_accum_exp (idm, src(PairHMM::IDM) + yTrans.lpTrans);
	    log_accum_exp (imi, src(PairHMM::IMI) + yTrans.lpTrans);
	  }
	}

	if (!xState.isNull() && !yState.isNull()) {
	  // xy-absorbing transitions into IMM
	  for (auto xt : xState.in) {
	    const ProfileTransition& xTrans = x.trans[xt];
	    for (auto yt : yState.in) {
	      const ProfileTransition& yTrans = y.trans[yt];
	      const XYCell& src = xyCell(xTrans.src,yTrans.src);

	      log_accum_exp (imm,
			     log_sum_exp (src(PairHMM::IMM) + hmm.imm_imm,
					  src(PairHMM::IMD) + hmm.imd_imm,
					  src(PairHMM::IDM) + hmm.idm_imm,
					  src(PairHMM::IMI) + hmm.imi_imm,
					  src(PairHMM::IIW) + hmm.iiw_imm)
			     + xTrans.lpTrans
			     + yTrans.lpTrans);
	    }
	  }

	  imm += computeLogProbAbsorb(i,j);

	} else {
	  if (xState.isNull() && i > 0) {
	    // x-nonabsorbing transitions in IMM
	    for (auto xt : xState.in) {
	      const ProfileTransition& xTrans = x.trans[xt];
	      log_accum_exp (imm, cell(xTrans.src,j,PairHMM::IMM) + xTrans.lpTrans);
	    }

	  } else if (yState.isNull() && j > 0) {
	    // y-nonabsorbing transitions in IMM
	    for (auto yt : yState.in) {
	      const ProfileTransition& yTrans = y.trans[yt];
	      log_accum_exp (imm, cell(i,yTrans.src,PairHMM::IMM) + yTrans.lpTrans);
	    }
	  }
	}
      }
    }
  }

  // transitions into EEE
  lpEnd = -numeric_limits<double>::infinity();
  for (auto xt : x.end().in) {
    const ProfileTransition& xTrans = x.trans[xt];
    for (auto yt : y.end().in) {
      const ProfileTransition& yTrans = y.trans[yt];
      log_accum_exp (lpEnd,
		     log_sum_exp (cell(xTrans.src,yTrans.src,PairHMM::IMM) + hmm.imm_eee,
				  cell(xTrans.src,yTrans.src,PairHMM::IMD) + hmm.imd_eee,
				  cell(xTrans.src,yTrans.src,PairHMM::IDM) + hmm.idm_eee,
				  cell(xTrans.src,yTrans.src,PairHMM::IMI) + hmm.imi_eee,
				  cell(xTrans.src,yTrans.src,PairHMM::IIW) + hmm.iiw_eee)
		     + xTrans.lpTrans
		     + yTrans.lpTrans);
    }
  }

  LogThisAt(6,"Forward log-likelihood is " << lpEnd << endl);
}

DPMatrix::CellCoords DPMatrix::sampleCell (const map<CellCoords,LogProb>& cellLogProb, random_engine& generator) const {
  double ptot = 0, lpmax = -numeric_limits<double>::infinity();
  for (auto& iter : cellLogProb)
    lpmax = max (lpmax, iter.second);
  for (auto& iter : cellLogProb)
    ptot += exp (iter.second - lpmax);
  uniform_real_distribution<double> dist (0, ptot);
  const double p0 = dist (generator);
  double p = p0;
  for (auto& iter : cellLogProb)
    if ((p -= exp (iter.second - lpmax)) <= 0) {
      LogThisAt(8,"Sampled " << cellName(iter.first) << " with probability " << exp(iter.second - lpmax) << "/" << ptot << endl);
      return iter.first;
    }
  for (auto& iter : cellLogProb)
    cerr << "Log P" << cellName(iter.first) << " = " << iter.second << endl;
  Abort ("%s fail (ptot=%g, p=%g)",__func__,ptot,p0);
  return CellCoords();
}

DPMatrix::CellCoords DPMatrix::bestCell (const map<CellCoords,LogProb>& cellLogProb) {
  CellCoords best;
  double pBest = -numeric_limits<double>::infinity();
  Assert (!cellLogProb.empty(), "%s traceback failure", __func__);
  for (auto& iter : cellLogProb)
    if (iter.second > pBest) {
      pBest = iter.second;
      best = iter.first;
    }
  return best;
}

ForwardMatrix::Path ForwardMatrix::sampleTrace (random_engine& generator) {
  Assert (lpEnd > -numeric_limits<double>::infinity(), "Forward likelihood is zero; traceback fail");
  
  Path path;
  path.push_back (endCell);

  map<CellCoords,LogProb> clp = sourceCells (endCell);
  CellCoords current;
  while (true) {
    current = sampleCell (clp, generator);
    LogThisAt(6,__func__ << " traceback at " << cellName(current) << " score " << cell(current) << endl);

    path.push_front (current);
    if (current.xpos == 0 && current.ypos == 0)
      break;
    clp = sourceCells (current);
  }

  return path;
}

ForwardMatrix::Path ForwardMatrix::bestTrace() {
  Assert (lpEnd > -numeric_limits<double>::infinity(), "Forward likelihood is zero; traceback fail");
  return bestTrace (endCell);
}

ForwardMatrix::Path ForwardMatrix::bestTrace (const CellCoords& end) {
  Path path;
  path.push_back (end);

  if (end.xpos > 0 || end.ypos > 0) {
    map<CellCoords,LogProb> clp = sourceCells (end);
    CellCoords current;
    while (true) {
      current = bestCell (clp);
      LogThisAt(6,__func__ << " traceback at " << cellName(current) << " score " << cell(current) << endl);
      
      path.push_front (current);
      if (current.xpos == 0 && current.ypos == 0)
	break;
      clp = sourceCells (current);
    }
  }

  return path;
}

AlignPath ForwardMatrix::bestAlignPath() {
  Path trace = bestTrace();
  return traceAlignPath (trace);
}

map<DPMatrix::CellCoords,LogProb> ForwardMatrix::sourceCells (const CellCoords& destCell) {
  map<CellCoords,LogProb> sc = sourceTransitions (destCell);
  for (auto& c_lp : sc)
    c_lp.second += cell (c_lp.first);
  return sc;
}

map<DPMatrix::CellCoords,LogProb> ForwardMatrix::sourceTransitions (const CellCoords& destCell) {
  map<CellCoords,LogProb> clp;
  const ProfileState& xState = x.state[destCell.xpos];
  const ProfileState& yState = y.state[destCell.ypos];

  switch (destCell.state) {
  case PairHMM::IMD:
  case PairHMM::IIW:
    if (xState.isNull()) {
      // x-nonabsorbing transitions in IMD, IIW
      if (destCell.xpos < xSize - 1)
	for (auto xt : xState.in)
	  clp[CellCoords(x.trans[xt].src,destCell.ypos,destCell.state)] = x.trans[xt].lpTrans;
    } else
      // x-absorbing transitions into IMD, IIW
      for (auto xt : xState.in)
	for (auto s : hmm.sources (destCell.state))
	  clp[CellCoords(x.trans[xt].src,destCell.ypos,s)] = hmm.lpTrans(s,destCell.state) + x.trans[xt].lpTrans;
    break;

  case PairHMM::IDM:
  case PairHMM::IMI:
    if (yState.isNull()) {
      // y-nonabsorbing transitions in IDM, IMI
      if (destCell.ypos < ySize - 1)
	for (auto yt : yState.in)
	  clp[CellCoords(destCell.xpos,y.trans[yt].src,destCell.state)] = y.trans[yt].lpTrans;
    } else
      // y-absorbing transitions into IDM, IMI
      for (auto yt : yState.in)
	for (auto s : hmm.sources (destCell.state))
	  clp[CellCoords(destCell.xpos,y.trans[yt].src,s)] = hmm.lpTrans(s,destCell.state) + y.trans[yt].lpTrans;
    break;

  case PairHMM::IMM:
    if (xState.isNull() && destCell.xpos > 0) {
      // x-nonabsorbing transitions in IMM
      if (destCell.xpos < xSize - 1)
	for (auto xt : xState.in)
	  clp[CellCoords(x.trans[xt].src,destCell.ypos,destCell.state)] = x.trans[xt].lpTrans;
    } else if (yState.isNull() && destCell.ypos > 0) {
      // y-nonabsorbing transitions in IMM
      if (destCell.ypos < ySize - 1)
	for (auto yt : yState.in)
	  clp[CellCoords(destCell.xpos,y.trans[yt].src,destCell.state)] = y.trans[yt].lpTrans;
    } else if (!xState.isNull() && !yState.isNull())
    // xy-absorbing transitions into IMM
      for (auto xt : xState.in)
	for (auto yt : yState.in)
	  for (auto s : hmm.sources (destCell.state))
	    clp[CellCoords(x.trans[xt].src,y.trans[yt].src,s)] = hmm.lpTrans(s,destCell.state) + x.trans[xt].lpTrans + y.trans[yt].lpTrans;
    break;

    // null transitions into EEE
  case PairHMM::EEE:
    if (destCell.xpos == xSize - 1 && destCell.ypos == ySize - 1)
      for (auto xt : x.end().in)
	for (auto yt : y.end().in)
	  for (auto s : hmm.sources (destCell.state))
	    clp[CellCoords(x.trans[xt].src,y.trans[yt].src,s)] = hmm.lpTrans (s,destCell.state) + x.trans[xt].lpTrans + y.trans[yt].lpTrans;
    break;

  default:
    Abort ("%s fail",__func__);
    break;
  }

  const LogProb lpAbs = lpCellEmitOrAbsorb (destCell);
  for (auto& src_lp : clp)
    src_lp.second += lpAbs;

  return clp;
}

DPMatrix::random_engine DPMatrix::newRNG() {
  return random_engine();
}

LogProb DPMatrix::lpCellEmitOrAbsorb (const CellCoords& c) {
  LogProb lp = 0;

  const ProfileState& xState = x.state[c.xpos];
  const ProfileState& yState = y.state[c.ypos];

  switch (c.state) {
  case PairHMM::IMD:
    if (!xState.isNull())
      lp = rootsubx[c.xpos];
    break;

  case PairHMM::IIW:
    if (!xState.isNull())
      lp = insx[c.xpos];
    break;

  case PairHMM::IDM:
    if (!yState.isNull())
      lp = rootsuby[c.ypos];
    break;

  case PairHMM::IMI:
    if (!yState.isNull())
      lp = insy[c.ypos];
    break;

  case PairHMM::IMM:
    if (!xState.isNull() && !yState.isNull())
      lp = computeLogProbAbsorb(c.xpos,c.ypos);

  default:
    break;
  }

  return lp;
}

string DPMatrix::cellName (const CellCoords& c) const {
  return string("(") + hmm.stateName(c.state,c.xpos==0,c.ypos==0) + "," + x.state[c.xpos].name + "," + y.state[c.ypos].name + ")";
}

bool DPMatrix::isAbsorbing (const CellCoords& c) const {
  return (c.state == PairHMM::IMM && !x.state[c.xpos].isNull() && !y.state[c.ypos].isNull())
    || (c.state == PairHMM::IMD && !x.state[c.xpos].isNull())
    || (c.state == PairHMM::IDM && !y.state[c.ypos].isNull());
}

bool DPMatrix::changesX (const CellCoords& c) const {
  return (c.state == PairHMM::IMM && !y.state[c.ypos].isNull())
    || c.state == PairHMM::IMD
    || c.state == PairHMM::IIW
    || c.state == PairHMM::EEE;
}

bool DPMatrix::changesY (const CellCoords& c) const {
  return (c.state == PairHMM::IMM && (!x.state[c.xpos].isNull() || y.state[c.ypos].isNull()))
    || c.state == PairHMM::IDM
    || c.state == PairHMM::IMI
    || c.state == PairHMM::EEE;
}

LogProb ForwardMatrix::eliminatedLogProbInsert (const CellCoords& cell) const {
  switch (cell.state) {
  case PairHMM::IIW:
    return x.state[cell.xpos].isNull() ? 0 : insx[cell.xpos];
    break;

  case PairHMM::IMI:
    return y.state[cell.ypos].isNull() ? 0 : insy[cell.ypos];
    break;

  case PairHMM::IMM:
  case PairHMM::IMD:
  case PairHMM::IDM:
  case PairHMM::EEE:
    return 0;
    break;
    
  default:
    Abort ("%s fail",__func__);
    break;
  }

  return -numeric_limits<double>::infinity();
}

ForwardMatrix::EffectiveTransition::EffectiveTransition()
  : lpPath (-numeric_limits<double>::infinity()),
    lpBestAlignPath (-numeric_limits<double>::infinity())
{ }

map<AlignRowIndex,SeqIdx> ForwardMatrix::cellSeqCoords (const CellCoords& c) const {
  map<AlignRowIndex,SeqIdx> coords = x.state[c.xpos].seqCoords;
  for (const auto& s_c : y.state[c.ypos].seqCoords)
    coords[s_c.first] = s_c.second;
  return coords;
}

AlignPath ForwardMatrix::cellAlignPath (const CellCoords& c) const {
  AlignPath alignPath;
  switch (c.state) {
  case PairHMM::IMM:
    if (!x.state[c.xpos].isNull() && !y.state[c.ypos].isNull())
      alignPath = alignPathUnion (x.state[c.xpos].alignPath, y.state[c.ypos].alignPath);
    else if (!x.state[c.xpos].isNull() && y.state[c.ypos].isNull())
      alignPath = y.state[c.ypos].alignPath;
    else // x.state[c.xpos].isNull()
      alignPath = x.state[c.xpos].alignPath;
    break;
  case PairHMM::IMD:
  case PairHMM::IIW:
    alignPath = x.state[c.xpos].alignPath;
    break;
  case PairHMM::IDM:
  case PairHMM::IMI:
    alignPath = y.state[c.ypos].alignPath;
    break;
  case PairHMM::EEE:
    break;
  default:
    Abort ("%s fail",__func__);
    break;
  }
  if (isAbsorbing(c))
    alignPath[parentRowIndex].push_back (true);
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

EigenCounts ForwardMatrix::transitionEigenCounts (const CellCoords& src, const CellCoords& dest) const {
  EigenCounts c;
  if (src.xpos != dest.xpos)
    c += x.getTrans(src.xpos,dest.xpos)->counts;
  if (src.ypos != dest.ypos)
    c += y.getTrans(src.ypos,dest.ypos)->counts;
  const bool xNull = x.state[dest.xpos].isNull();
  const bool yNull = y.state[dest.ypos].isNull();
  switch (dest.state) {
  case PairHMM::IMM:
    if (!xNull && !yNull)
      c.indelCounts.matchTime += hmm.l.t + hmm.r.t;
    break;
  case PairHMM::IMD:
    if (!xNull) {
      c.indelCounts.matchTime += hmm.l.t;
      if (src.state == dest.state)
	c.indelCounts.delExt += 1;
      else {
	c.indelCounts.del += 1;
	c.indelCounts.delTime += hmm.r.t;
      }
    }
    break;
  case PairHMM::IIW:
    if (!xNull) {
      if (src.state == dest.state)
	c.indelCounts.insExt += 1;
      else
	c.indelCounts.ins += 1;
    }
    break;
  case PairHMM::IDM:
    if (!yNull) {
      c.indelCounts.matchTime += hmm.r.t;
      if (src.state == dest.state)
	c.indelCounts.delExt += 1;
      else {
	c.indelCounts.del += 1;
	c.indelCounts.delTime += hmm.l.t;
      }
    }
    break;
  case PairHMM::IMI:
    if (!yNull) {
      if (src.state == dest.state)
	c.indelCounts.insExt += 1;
      else
	c.indelCounts.ins += 1;
    }
    break;
  default:
    break;
  }
  return c;
}

AlignPath ForwardMatrix::traceAlignPath (const Path& path) const {
  AlignPath p;
  const vector<CellCoords> pv (path.begin(), path.end());
  for (size_t n = 0; n < pv.size() - 1; ++n)
    p = alignPathConcat (p, cellAlignPath(pv[n]), transitionAlignPath(pv[n],pv[n+1]));
  p = alignPathConcat (p, cellAlignPath(pv.back()));

  const AlignRowIndex rows = p.size();
  const AlignColIndex cols = alignPathColumns (p);  // this will also test if alignment is flush
  LogThisAt(2,"Converted Forward matrix trace into alignment with " << rows << " rows and " << cols << " columns" << endl);

  return p;
}

Profile ForwardMatrix::makeProfile (const set<CellCoords>& cells, ProfilingStrategy strategy) {
  Profile prof (alphSize);
  prof.name = Tree::pairParentName (x.name, hmm.l.t, y.name, hmm.r.t);
  prof.meta["node"] = to_string(parentRowIndex);

  Assert (cells.find (startCell) != cells.end(), "Missing SSS");
  Assert (cells.find (endCell) != cells.end(), "Missing EEE");

  // build states
  // retain only start, end, and absorbing cells
  map<CellCoords,ProfileStateIndex> profStateIndex;
  map<CellCoords,int> outgoingTransitionCount;

  for (const auto& dest : cells)
    for (const auto& src_lp : sourceTransitions (dest))
      ++outgoingTransitionCount[src_lp.first];

  for (const auto& c : cells)
    if (isAbsorbing(c)
	|| c == startCell
	|| c == endCell
	|| ((strategy & CollapseChains) != 0 && outgoingTransitionCount[c] > 1)
	|| (strategy & CollapseChains) == 0) {
      // cell is to be retained
      profStateIndex[c] = prof.state.size();
      prof.state.push_back (ProfileState());
      if (isAbsorbing(c))
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
	  break;
	}
      prof.state.back().alignPath = cellAlignPath(c);
      prof.state.back().seqCoords = cellSeqCoords(c);
      prof.state.back().name = cellName(c);
      prof.state.back().meta["fwdLogProb"] = to_string(c.state == PairHMM::EEE ? lpEnd : cell(c.xpos,c.ypos,c.state));
    }

  // Calculate log-probabilities of "effective transitions" from cells to retained-cells.
  // Each effective transition represents a sum over paths.
  // A path is either a single direct transition (from source cell to destination retained-cell),
  // or a series of transitions starting from the source cell,
  // passing through one or more eliminated-cells, and stopping at the destination retained-cell.
  map<CellCoords,map<ProfileStateIndex,EffectiveTransition> > effTrans;  // effTrans[srcCell][destStateIdx]
  for (auto iter = cells.crbegin(); iter != cells.crend(); ++iter) {
    const CellCoords& cell = *iter;
    const map<CellCoords,LogProb>& slp = sourceTransitions (cell);
    const LogProb cellLogProbInsert = eliminatedLogProbInsert (cell);
    if (profStateIndex.find(cell) != profStateIndex.end()) {
      // cell is to be retained. Incoming & outgoing paths can be kept separate
      const ProfileStateIndex cellIdx = profStateIndex[cell];
      for (const auto slpIter : slp) {
	const CellCoords& src = slpIter.first;
	const LogProb srcCellLogProbTrans = slpIter.second;
	EffectiveTransition& eff = effTrans[src][cellIdx];
	eff.lpPath = eff.lpBestAlignPath = srcCellLogProbTrans + cellLogProbInsert;
	eff.bestAlignPath = transitionAlignPath(src,cell);
	if (strategy & CountEvents)
	  eff.counts = transitionEigenCounts(src,cell);
      }
    } else {
      // cell is to be eliminated. Connect incoming transitions & outgoing paths, summing cell out
      const auto& cellEffTrans = effTrans[cell];
      const AlignPath& cap = cellAlignPath (cell);
      EigenCounts cellCounts, srcCellCounts;
      if ((strategy & CountEvents) != 0 && sumProd != NULL)
	cellCounts = cachedCellEigenCounts (cell, *sumProd);
      for (auto slpIter : slp) {
	const CellCoords& src = slpIter.first;
	const LogProb srcCellLogProbTrans = slpIter.second;
	if (strategy & CountEvents)
	  srcCellCounts = transitionEigenCounts (src, cell) + cellCounts;
	auto& srcEffTrans = effTrans[src];
	for (const auto cellEffTransIter : cellEffTrans) {
	  const ProfileStateIndex& destIdx = cellEffTransIter.first;
	  const EffectiveTransition& cellDestEffTrans = cellEffTransIter.second;
	  EffectiveTransition& srcDestEffTrans = srcEffTrans[destIdx];
	  const LogProb lpPath = srcCellLogProbTrans + cellLogProbInsert + cellDestEffTrans.lpPath;
	  log_accum_exp (srcDestEffTrans.lpPath, lpPath);
	  if (strategy & CountEvents) {
	    const double ppPath = exp (lpPath - srcDestEffTrans.lpPath);
	    srcDestEffTrans.counts *= 1 - ppPath;
	    srcDestEffTrans.counts += (srcCellCounts + cellDestEffTrans.counts) * ppPath;
	  }
	  // we also want to keep the best single alignment consistent w/this set of paths
	  const LogProb srcDestLogProbBestAlignPath = srcCellLogProbTrans + cellLogProbInsert + cellDestEffTrans.lpBestAlignPath;
	  if (srcDestLogProbBestAlignPath > srcDestEffTrans.lpBestAlignPath) {
	    srcDestEffTrans.lpBestAlignPath = srcDestLogProbBestAlignPath;
	    srcDestEffTrans.bestAlignPath = alignPathConcat (transitionAlignPath(src,cell), cap, cellDestEffTrans.bestAlignPath);
	  }
	}
      }
    }
  }

  // populate outgoing & incoming transitions for each state
  for (const auto profStateIter : profStateIndex) {
    const CellCoords& cell = profStateIter.first;
    const ProfileStateIndex srcIdx = profStateIter.second;
    vguard<ProfileTransitionIndex>& srcNullOut = prof.state[srcIdx].nullOut;
    vguard<ProfileTransitionIndex>& srcAbsorbOut = prof.state[srcIdx].absorbOut;
    for (const auto effTransIter : effTrans[cell]) {
      const ProfileStateIndex destIdx = effTransIter.first;
      const EffectiveTransition& srcDestEffTrans = effTransIter.second;
      vguard<ProfileTransitionIndex>& destIn = prof.state[destIdx].in;
      const ProfileTransitionIndex transIdx = prof.trans.size();
      ProfileTransition trans;
      trans.src = srcIdx;
      trans.dest = destIdx;
      trans.lpTrans = srcDestEffTrans.lpPath;
      trans.alignPath = srcDestEffTrans.bestAlignPath;
      if (strategy & CountEvents)
	trans.counts = srcDestEffTrans.counts;
      prof.trans.push_back (trans);
      (prof.state[destIdx].isNull() ? srcNullOut : srcAbsorbOut).push_back (transIdx);
      destIn.push_back (transIdx);
    }
  }

  prof.seq = x.seq;
  prof.seq.insert (y.seq.begin(), y.seq.end());
  
  return prof;
}

Profile ForwardMatrix::sampleProfile (random_engine& generator, size_t profileSamples, size_t maxCells, ProfilingStrategy strategy) {
  map<CellCoords,size_t> cellCount;

  Require ((strategy & IncludeBestTrace) || profileSamples > 0, "Must allow at least one sample path in the profile");

  size_t nTraces = 0;
  if (strategy & IncludeBestTrace) {
    const Path best = bestTrace();
    for (auto& c : best)
      cellCount[c] = 2;  // avoid dropping these cells
    ++nTraces;
  }
  for (size_t n = 0; n < profileSamples && (maxCells == 0 || cellCount.size() < maxCells); ++n) {
    const Path sampled = sampleTrace (generator);
    if (LoggingThisAt(5)) {
      LogThisAt(5,"Trace #" << n+1 << ":");
      for (auto& c : sampled)
	LogThisAt(5," " << cellName(c));
      LogThisAt(5,endl);
    }
    for (auto& c : sampled)
      ++cellCount[c];
    ++nTraces;
  }
  set<CellCoords> profCells;
  const size_t threshold = (nTraces > 1 && maxCells > 0 && cellCount.size() >= maxCells) ? 2 : 1;
  for (const auto& cc : cellCount)
    if (cc.second >= threshold)
      profCells.insert (cc.first);
  return makeProfile (profCells, strategy);
}

Profile ForwardMatrix::bestProfile (ProfilingStrategy strategy) {
  const Path best = bestTrace();
  const set<CellCoords> profCells (best.begin(), best.end());
  return makeProfile (profCells, strategy);
}

EigenCounts ForwardMatrix::cellEigenCounts (const CellCoords& cell, SumProduct& sumProd) const {
  EigenCounts c (hmm.alphabetSize());
  accumulateEigenCounts (c, cell, sumProd);
  return c;
}

EigenCounts ForwardMatrix::cachedCellEigenCounts (const CellCoords& cell, SumProduct& sumProd) {
  if (!isAbsorbing (cell)) {
    if (changesX (cell)) {
      if (xInsertCounts.find (cell.xpos) == xInsertCounts.end())
	xInsertCounts[cell.xpos] = cellEigenCounts (cell, sumProd);
      return xInsertCounts[cell.xpos];

    } else if (changesY (cell)) {
      if (yInsertCounts.find (cell.ypos) == yInsertCounts.end())
	yInsertCounts[cell.ypos] = cellEigenCounts (cell, sumProd);
      return yInsertCounts[cell.ypos];
    }
  }

  return cellEigenCounts (cell, sumProd);
}

void ForwardMatrix::accumulateEigenCounts (EigenCounts& counts, const CellCoords& cell, SumProduct& sumProd, double weight) const {
  LogThisAt(9,"Accumulating event counts for cell " << cellName(cell) << endl);
  const auto col = getAlignmentColumn (cell);
  if (col.size()) {
    sumProd.initColumn (col);
    sumProd.fillUp();
    sumProd.fillDown();
    sumProd.accumulateEigenCounts (counts.rootCount, counts.eigenCount, weight);
  }
}

void ForwardMatrix::accumulateCachedEigenCounts (EigenCounts& counts, const CellCoords& cell, SumProduct& sumProd, double weight) {
  if (!isAbsorbing(cell) && (changesX(cell) || changesY(cell)))
    counts += cachedCellEigenCounts (cell, sumProd) * weight;
  else
    accumulateEigenCounts (counts, cell, sumProd, weight);
}

map<AlignRowIndex,char> ForwardMatrix::getAlignmentColumn (const CellCoords& cell) const {
  map<AlignRowIndex,char> col;
  if (cell.xpos > 0 && cell.ypos > 0 && cell.xpos < xSize - 1 && cell.ypos < ySize - 1)
    switch (cell.state) {
    case PairHMM::IMM:
      if (!x.state[cell.xpos].isNull() && !y.state[cell.ypos].isNull()) {
	col = x.alignColumn (cell.xpos);
	const auto yCol = y.alignColumn (cell.ypos);
	col.insert (yCol.begin(), yCol.end());
	col[parentRowIndex] = Alignment::wildcardChar;
      }
      break;
    case PairHMM::IMD:
      if (!x.state[cell.xpos].isNull()) {
	col = x.alignColumn (cell.xpos);
	col[parentRowIndex] = Alignment::wildcardChar;
      }
      break;
    case PairHMM::IDM:
      if (!y.state[cell.ypos].isNull()) {
	col = y.alignColumn (cell.ypos);
	col[parentRowIndex] = Alignment::wildcardChar;
      }
      break;
    case PairHMM::IIW:
      if (!x.state[cell.xpos].isNull())
	col = x.alignColumn (cell.xpos);
      break;
    case PairHMM::IMI:
      if (!y.state[cell.ypos].isNull())
	col = y.alignColumn (cell.ypos);
      break;
    default:
      break;
    }
  return col;
}

BackwardMatrix::BackwardMatrix (ForwardMatrix& fwd, double minPostProb)
  : DPMatrix (fwd.x, fwd.y, fwd.hmm, fwd.envelope),
    fwd (fwd)
{
  lpEnd = 0;
  // transitions into EEE
  for (auto xt : x.end().in) {
    const ProfileTransition& xTrans = x.trans[xt];
    for (auto yt : y.end().in) {
      const ProfileTransition& yTrans = y.trans[yt];
      if (envelope.inRange (xClosestLeafPos[xTrans.src], yClosestLeafPos[yTrans.src])) {
	XYCell& src = xyCell(xTrans.src,yTrans.src);
      
	src(PairHMM::IMM) = xTrans.lpTrans + yTrans.lpTrans + hmm.imm_eee;
	src(PairHMM::IMD) = xTrans.lpTrans + yTrans.lpTrans + hmm.imd_eee;
	src(PairHMM::IDM) = xTrans.lpTrans + yTrans.lpTrans + hmm.idm_eee;
	src(PairHMM::IMI) = xTrans.lpTrans + yTrans.lpTrans + hmm.imi_eee;
	src(PairHMM::IIW) = xTrans.lpTrans + yTrans.lpTrans + hmm.iiw_eee;
      }
    }
  }

  const LogProb fwdEnd = fwd.lpEnd;
  const LogProb lppThreshold = log(minPostProb);
  const auto states = hmm.states();

  ProgressLog (plog, 4);
  plog.initProgress ("Backward algorithm (%s vs %s)", x.name.c_str(), y.name.c_str());

  for (int i = xSize - 2; i >= 0; --i) {
    const ProfileState& xState = x.state[i];

    plog.logProgress ((xSize - 2 - i) / (double) (xSize - 2), "state %d/%d", xSize - 1 - i, xSize);

    for (int j = ySize - 2; j >= 0; --j) {
      const ProfileState& yState = y.state[j];

      if (envelope.inRange (xClosestLeafPos[i], yClosestLeafPos[j])) {

	XYCell& src = xyCell(i,j);
	double& imm = src(PairHMM::IMM);
	double& imd = src(PairHMM::IMD);
	double& idm = src(PairHMM::IDM);
	double& imi = src(PairHMM::IMI);
	double& iiw = src(PairHMM::IIW);

	// xy-absorbing transitions into IMM
	for (auto xt : xState.absorbOut) {
	  const ProfileTransition& xTrans = x.trans[xt];
	  for (auto yt : yState.absorbOut) {
	    const ProfileTransition& yTrans = y.trans[yt];
	    const LogProb dest_imm = xTrans.lpTrans + yTrans.lpTrans + computeLogProbAbsorb(xTrans.dest,yTrans.dest) + cell(xTrans.dest,yTrans.dest,PairHMM::IMM);
	    
	    log_accum_exp (imm, hmm.imm_imm + dest_imm);
	    log_accum_exp (imd, hmm.imd_imm + dest_imm);
	    log_accum_exp (idm, hmm.idm_imm + dest_imm);
	    log_accum_exp (imi, hmm.imi_imm + dest_imm);
	    log_accum_exp (iiw, hmm.iiw_imm + dest_imm);
	  }
	}

	// x-absorbing transitions into IMD, IIW
	for (auto xt : xState.absorbOut) {
	  const ProfileTransition& xTrans = x.trans[xt];
	  const XYCell& dest = xyCell(xTrans.dest,j);
	  const LogProb dest_imd = xTrans.lpTrans + rootsubx[xTrans.dest] + dest(PairHMM::IMD);
	  const LogProb dest_iiw = xTrans.lpTrans + insx[xTrans.dest] + dest(PairHMM::IIW);

	  log_accum_exp (imm, hmm.imm_imd + dest_imd);
	  log_accum_exp (imd, hmm.imd_imd + dest_imd);
	  log_accum_exp (idm, hmm.idm_imd + dest_imd);
	  log_accum_exp (imi, hmm.imi_imd + dest_imd);

	  log_accum_exp (imm, hmm.imm_iiw + dest_iiw);
	  log_accum_exp (imi, hmm.imi_iiw + dest_iiw);
	  log_accum_exp (iiw, hmm.iiw_iiw + dest_iiw);
	}

	// y-absorbing transitions into IDM, IMI
	for (auto yt : yState.absorbOut) {
	  const ProfileTransition& yTrans = y.trans[yt];
	  const XYCell& dest = xyCell(i,yTrans.dest);
	  const LogProb dest_idm = yTrans.lpTrans + rootsuby[yTrans.dest] + dest(PairHMM::IDM);
	  const LogProb dest_imi = yTrans.lpTrans + insy[yTrans.dest] + dest(PairHMM::IMI);

	  log_accum_exp (imm, hmm.imm_idm + dest_idm);
	  log_accum_exp (imd, hmm.imd_idm + dest_idm);
	  log_accum_exp (idm, hmm.idm_idm + dest_idm);
	  log_accum_exp (iiw, hmm.iiw_idm + dest_idm);

	  log_accum_exp (imm, hmm.imm_imi + dest_imi);
	  log_accum_exp (imi, hmm.imi_imi + dest_imi);
	}

	// x-nonabsorbing transitions in IMD, IIW, IMM
	if (!yState.isNull())
	  for (auto xt : xState.nullOut) {
	    const ProfileTransition& xTrans = x.trans[xt];
	    const XYCell& dest = xyCell(xTrans.dest,j);
	    log_accum_exp (imd, xTrans.lpTrans + dest(PairHMM::IMD));
	    log_accum_exp (iiw, xTrans.lpTrans + dest(PairHMM::IIW));
	    log_accum_exp (imm, xTrans.lpTrans + dest(PairHMM::IMM));
	  }
	
	// y-nonabsorbing transitions in IDM, IMI, IMM
	for (auto yt : yState.nullOut) {
	  const ProfileTransition& yTrans = y.trans[yt];
	  const XYCell& dest = xyCell(i,yTrans.dest);
	  log_accum_exp (idm, yTrans.lpTrans + dest(PairHMM::IDM));
	  log_accum_exp (imi, yTrans.lpTrans + dest(PairHMM::IMI));
	  log_accum_exp (imm, yTrans.lpTrans + dest(PairHMM::IMM));
	}

	// track cells above posterior probability threshold
	const XYCell& fwdSrc = fwd.xyCell(i,j);
	for (auto s : states) {
	  const LogProb lpp = src(s) + fwdSrc(s) - fwdEnd;
	  if (lpp >= lppThreshold)
	    bestCells.push (CellPostProb (i, j, s, lpp));
	}
      }
    }
  }

  LogThisAt(6,"Backward log-likelihood is " << lpStart() << endl);
  if (gsl_fcmp (lpStart(), fwd.lpEnd, FWD_BACK_ERROR_TOLERANCE) != 0) {
    fwd.slowFillTest();
    slowFillTest();
    sourceDestTransTest();
    Abort ("Forward log-likelihood is %g, Backward log-likelihood is %g", fwd.lpEnd, lpStart());
  }
}

void ForwardMatrix::slowFillTest() {
  auto states = hmm.states();
  states.push_back (PairHMM::EEE);
  size_t nCells = 0, nTrans = 0;
  for (int i = 0; i < xSize; ++i)
    for (int j = 0; j < ySize; ++j)
      if (envelope.inRange (xClosestLeafPos[i], yClosestLeafPos[j]))
	for (auto s : states) {
	  const bool atStart = s == PairHMM::SSS && i == 0 && j == 0;
	  const bool atEnd = s == PairHMM::EEE && i == xSize-1 && j == ySize-1;
	  if ((i < xSize-1 && j < ySize-1 && s != PairHMM::EEE) || atEnd) {
	    ++nCells;
	    const CellCoords destCell (i, j, s);
	    const LogProb lpDestCell = atEnd ? lpEnd : cell(destCell);
	    LogProb lp = atStart ? 0 : -numeric_limits<double>::infinity();
	    for (auto src_lp : sourceTransitions (destCell))
	      if (src_lp.second > -numeric_limits<double>::infinity()) {
		log_accum_exp_slow (lp, src_lp.second + cell(src_lp.first));
		++nTrans;
	      }
	    Test (gsl_fcmp (lp, lpDestCell, FWD_BACK_ERROR_TOLERANCE) == 0, "Forward cell %s score (%g) doesn't match slow computation (%g)", cellName(destCell).c_str(), lpDestCell, lp);
	  }
	}
  LogThisAt(6,"Forward slow fill test: iterated over " << nCells << " cells and " << nTrans << " transitions" << endl);
}

void BackwardMatrix::slowFillTest() {
  const auto states = hmm.states();
  size_t nCells = 0, nTrans = 0;
  for (int i = xSize - 2; i >= 0; --i)
    for (int j = ySize - 2; j >= 0; --j)
      if (envelope.inRange (xClosestLeafPos[i], yClosestLeafPos[j]))
	for (auto s : states) {
	  ++nCells;
	  const CellCoords srcCell (i, j, s);
	  LogProb lp = -numeric_limits<double>::infinity();
	  for (auto dest_lp : destTransitions (srcCell))
	    if (dest_lp.second > -numeric_limits<double>::infinity()) {
	      log_accum_exp_slow (lp, dest_lp.second + (dest_lp.first.state == PairHMM::EEE ? 0 : cell(dest_lp.first)));
	      ++nTrans;
	    }
	  Test (gsl_fcmp (lp, cell(srcCell), FWD_BACK_ERROR_TOLERANCE) == 0, "Backward cell %s score (%g) doesn't match slow computation (%g)", cellName(srcCell).c_str(), cell(srcCell), lp);
	}
  ++nCells;  // account for endCell
  LogThisAt(6,"Backward slow fill test: iterated over " << nCells << " cells and " << nTrans << " transitions" << endl);
}

void BackwardMatrix::sourceDestTransTest() {
  const auto states = hmm.states();
  for (int i = 0; i < xSize; ++i)
    for (int j = 0; j < ySize; ++j)
      if (envelope.inRange (xClosestLeafPos[i], yClosestLeafPos[j]))
	for (auto s : states) {
	  const CellCoords cell (i, j, s);
	  for (const auto& src_lp : fwd.sourceTransitions (cell))
	    if (src_lp.second > -numeric_limits<double>::infinity()) {
	      auto destTrans = destTransitions (src_lp.first);
	      Test (gsl_fcmp (src_lp.second, destTrans[cell], FWD_BACK_ERROR_TOLERANCE) == 0, "Forward (%g) & Backward (%g) transitions between %s and %s don't match", src_lp.second, destTrans[cell], cellName(src_lp.first).c_str(), cellName(cell).c_str());
	    }
	  for (const auto& dest_lp : destTransitions (cell))
	    if (dest_lp.second > -numeric_limits<double>::infinity()) {
	      auto srcTrans = fwd.sourceTransitions (dest_lp.first);
	      Test (gsl_fcmp (dest_lp.second, srcTrans[cell], FWD_BACK_ERROR_TOLERANCE) == 0, "Forward (%g) & Backward (%g) transitions between %s and %s don't match", srcTrans[cell], dest_lp.second, cellName(cell).c_str(), cellName(dest_lp.first).c_str());
	    }
	}
}

double BackwardMatrix::cellPostProb (const CellCoords& c) const {
  return exp (fwd.cell(c) + cell(c) - fwd.lpEnd);
}

double BackwardMatrix::transPostProb (const CellCoords& src, const CellCoords& dest) const {
  const auto srcTrans = fwd.sourceTransitions (dest);
  if (srcTrans.find (src) != srcTrans.end())
    return exp (fwd.cell(src) + srcTrans.at(src) + cell(dest) - fwd.lpEnd);
  return 0;
}

EigenCounts BackwardMatrix::getCounts() const {
  EigenCounts counts (hmm.alphabetSize());
  const auto states = hmm.states();

  ProgressLog (plog, 4);
  plog.initProgress ("Forward-Backward counts (%s vs %s)", x.name.c_str(), y.name.c_str());

  for (ProfileStateIndex i = 0; i < xSize - 1; ++i) {
    const ProfileState& xState = x.state[i];
    plog.logProgress (i / (double) (xSize - 2), "state %d/%d", i + 1, xSize);

    for (ProfileStateIndex j = 0; j < ySize - 1; ++j) {
      const ProfileState& yState = y.state[j];

      if (envelope.inRange (xClosestLeafPos[i], yClosestLeafPos[j])) {
	for (auto s : states) {
	  const CellCoords dest (i, j, s);
	  const LogProb lpDest = cell(dest);
	  if (fwd.sumProd)
	    fwd.accumulateCachedEigenCounts (counts, dest, *fwd.sumProd, exp (fwd.cell(dest) + lpDest - fwd.lpEnd));
	  const auto srcTrans = fwd.sourceTransitions (dest);
	  for (auto& src_lp : srcTrans)
	    counts += fwd.transitionEigenCounts (src_lp.first, dest) * exp (fwd.cell(src_lp.first) + src_lp.second + lpDest - fwd.lpEnd);
	}
      }
    }
  }

  return counts;
}

map<DPMatrix::CellCoords,LogProb> BackwardMatrix::destCells (const CellCoords& srcCell) {
  map<CellCoords,LogProb> clp = destTransitions (srcCell);
  for (auto& c_lp : clp)
    if (c_lp.first.state != PairHMM::EEE)
      c_lp.second += cell (c_lp.first);
  return clp;
}

map<DPMatrix::CellCoords,LogProb> BackwardMatrix::destTransitions (const CellCoords& srcCell) {
  map<CellCoords,LogProb> clp;
  const ProfileState& xState = x.state[srcCell.xpos];
  const ProfileState& yState = y.state[srcCell.ypos];
  
  // xy-absorbing transitions into IMM
  for (auto xt : xState.absorbOut) {
    const ProfileTransition& xTrans = x.trans[xt];
    for (auto yt : yState.absorbOut) {
      const ProfileTransition& yTrans = y.trans[yt];
      clp[CellCoords(xTrans.dest,yTrans.dest,PairHMM::IMM)] = hmm.lpTrans(srcCell.state,PairHMM::IMM) + xTrans.lpTrans + yTrans.lpTrans;
    }
  }

  // x-absorbing transitions into IMD, IIW
  for (auto xt : xState.absorbOut) {
    const ProfileTransition& xTrans = x.trans[xt];
    clp[CellCoords(xTrans.dest,srcCell.ypos,PairHMM::IMD)] = hmm.lpTrans(srcCell.state,PairHMM::IMD) + xTrans.lpTrans;
    clp[CellCoords(xTrans.dest,srcCell.ypos,PairHMM::IIW)] = hmm.lpTrans(srcCell.state,PairHMM::IIW) + xTrans.lpTrans;
  }

  // y-absorbing transitions into IDM, IMI
  for (auto yt : yState.absorbOut) {
    const ProfileTransition& yTrans = y.trans[yt];
    clp[CellCoords(srcCell.xpos,yTrans.dest,PairHMM::IDM)] = hmm.lpTrans(srcCell.state,PairHMM::IDM) + yTrans.lpTrans;
    clp[CellCoords(srcCell.xpos,yTrans.dest,PairHMM::IMI)] = hmm.lpTrans(srcCell.state,PairHMM::IMI) + yTrans.lpTrans;
  }

  // x-nonabsorbing transitions in IMD, IIW, IMM
  if (!yState.isNull() && (srcCell.state == PairHMM::IMD || srcCell.state == PairHMM::IIW || srcCell.state == PairHMM::IMM))
    for (auto xt : xState.nullOut) {
      const ProfileTransition& xTrans = x.trans[xt];
      if (xTrans.dest != xSize - 1)
	clp[CellCoords(xTrans.dest,srcCell.ypos,srcCell.state)] = xTrans.lpTrans;
    }

  // y-nonabsorbing transitions in IDM, IMI, IMM
  if (srcCell.state == PairHMM::IDM || srcCell.state == PairHMM::IMI || srcCell.state == PairHMM::IMM)
    for (auto yt : yState.nullOut) {
      const ProfileTransition& yTrans = y.trans[yt];
      if (yTrans.dest != ySize - 1)
	clp[CellCoords(srcCell.xpos,yTrans.dest,srcCell.state)] = yTrans.lpTrans;
    }

  // add in transitions to EEE
  for (auto xt : xState.nullOut) {
    const ProfileTransition& xTrans = x.trans[xt];
    if (xTrans.dest == xSize - 1)
      for (auto yt : yState.nullOut) {
	const ProfileTransition& yTrans = y.trans[yt];
	if (yTrans.dest == ySize - 1)
	  clp[CellCoords(xTrans.dest,yTrans.dest,PairHMM::EEE)] = xTrans.lpTrans + yTrans.lpTrans + hmm.lpTrans(srcCell.state,PairHMM::EEE);
      }
  }

  for (auto& dest_lp : clp)
    dest_lp.second += lpCellEmitOrAbsorb (dest_lp.first);

  return clp;
}

BackwardMatrix::Path BackwardMatrix::bestTrace (const CellCoords& traceStart) {
  Path path;

  CellCoords current = traceStart;
  while (current.xpos < xSize - 1 && current.ypos < ySize - 1) {
    map<CellCoords,LogProb> clp = destCells (current);
    current = bestCell (clp);
    LogThisAt(6,__func__ << " traceforward at " << cellName(current) << " score " << cell(current) << endl);
    path.push_back (current);
  }

  path.push_back (endCell);
  return path;
}

Profile BackwardMatrix::buildProfile (size_t maxCells, ProfilingStrategy strategy) {
  set<CellCoords> cells;
  if (bestCells.empty() || (strategy & IncludeBestTrace)) {
    const list<CellCoords> fwdBestTrace = fwd.bestTrace();
    cells.insert (fwdBestTrace.begin(), fwdBestTrace.end());
  }
  while ((maxCells == 0 || cells.size() < maxCells) && !bestCells.empty()) {
    const CellCoords& best = bestCells.top();
    if (cells.count (best))
      bestCells.pop();
    else {
      LogThisAt(5,"Starting traceback/forward from " << cellName(best) << endl);
      const list<CellCoords> fwdTrace = fwd.bestTrace(best), backTrace = bestTrace(best);
      list<CellCoords> newCells;
      for (auto cellIter = fwdTrace.rbegin(); cellIter != fwdTrace.rend(); ++cellIter)
	if (cells.count (*cellIter))
	  break;
	else
	  newCells.push_back (*cellIter);
      for (const auto& cell : backTrace)
	if (cells.count (cell))
	  break;
	else
	  newCells.push_back (cell);
      if (maxCells > 0 && cells.size() > 0 && cells.size() + newCells.size() > maxCells)
	break;
      if (LoggingThisAt(6)) {
	LogThisAt(6,"Adding the following cells to profile:");
	for (const auto& c : newCells)
	  LogThisAt(6," " << cellName(c));
	LogThisAt(6,endl);
      }
      cells.insert (newCells.begin(), newCells.end());
    }
  }
  return fwd.makeProfile (cells, strategy);
}
