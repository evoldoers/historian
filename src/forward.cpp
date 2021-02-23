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
    xEmpty (x.isEmpty()),
    yEmpty (y.isEmpty()),
    xSize (x.size()),
    ySize (y.size()),
    subx (x.leftMultiply (hmm.l.subMat)),
    suby (y.leftMultiply (hmm.r.subMat)),
    cellStorage (x.size()),
    absorbScratch (hmm.components(), vguard<LogProb> (hmm.alphabetSize())),
    insx (x.size(), -numeric_limits<double>::infinity()),
    insy (y.size(), -numeric_limits<double>::infinity()),
    rootsubx (x.size(), -numeric_limits<double>::infinity()),
    rootsuby (y.size(), -numeric_limits<double>::infinity()),
    startCell (0, 0, PairHMM::SSS),
    endCell (xSize - 1, ySize - 1, PairHMM::EEE),
    envelope (env),
    xClosestLeafPos (xSize, 0),
    yClosestLeafPos (ySize, 0),
    xNearStart (xSize, false),
    yNearEnd (ySize, false)
{
  if (env.initialized()) {
    for (ProfileStateIndex i = 1; i < xSize; ++i)
      xClosestLeafPos[i] = x.state[i].seqCoords.at(env.row1);

    for (ProfileStateIndex j = 1; j < ySize; ++j)
      yClosestLeafPos[j] = y.state[j].seqCoords.at(env.row2);
  }

  for (ProfileStateIndex i = 1; i < xSize - 1; ++i)
    if (!x.state[i].isNull())
      for (int cpt = 0; cpt < components(); ++cpt) {
	log_accum_exp (insx[i], hmm.logl.logCptWeight[cpt] + logInnerProduct (hmm.logl.logInsProb[cpt], x.state[i].lpAbsorb[cpt]));
	log_accum_exp (rootsubx[i], logInnerProduct (hmm.logRoot[cpt], subx.state[i].lpAbsorb[cpt]));
      }

  for (ProfileStateIndex j = 1; j < ySize - 1; ++j)
    if (!y.state[j].isNull())
      for (int cpt = 0; cpt < components(); ++cpt) {
	log_accum_exp (insy[j], hmm.logr.logCptWeight[cpt] + logInnerProduct (hmm.logr.logInsProb[cpt], y.state[j].lpAbsorb[cpt]));
	log_accum_exp (rootsuby[j], logInnerProduct (hmm.logRoot[cpt], suby.state[j].lpAbsorb[cpt]));
      }

  xNearStart[0] = true;
  for (ProfileStateIndex i = 0; i < xSize; ++i)
    if (xNearStart[i])
      for (const auto& t: x.state[i].nullOut)
	xNearStart[x.trans[t].dest] = true;
  
  for (auto yt : y.end().in)
    yNearEnd[y.trans[yt].src] = true;
}

ForwardMatrix::ForwardMatrix (const Profile& x, const Profile& y, const PairHMM& hmm, AlignRowIndex parentRowIndex, const GuideAlignmentEnvelope& env, SumProduct* sumProd)
  : DPMatrix (x, y, hmm, env),
    parentRowIndex (parentRowIndex),
    sumProd (sumProd)
{
  lpStart() = 0;

  ProgressLog (plog, 5);
  plog.initProgress ("Forward algorithm (%s vs %s)", x.name.c_str(), y.name.c_str());

  for (ProfileStateIndex i = 0; i < xSize - 1; ++i) {
    const ProfileState& xState = x.state[i];

    plog.logProgress (i / (double) (xSize - 2), "state %d/%d", i + 1, xSize);

    for (ProfileStateIndex j = 0; j < ySize - 1; ++j) {
      const ProfileState& yState = y.state[j];

      if (inEnvelope(i,j)) {

	XYCell& dest = xyCell(i,j);
	double& imm = dest(PairHMM::IMM);
	double& imd = dest(PairHMM::IMD);
	double& idm = dest(PairHMM::IDM);
	double& imi = dest(PairHMM::IMI);
	double& iiw = dest(PairHMM::IIW);

	if (!xState.isNull()) {
	  // x-absorbing transitions into IMD, IIW
	  if (yState.isReady() || yEmpty) {
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
	  }

	} else {  // xState.isNull()
	  // x-nonabsorbing transitions in IMD, IIW
	  if (yState.isReady() || yEmpty)
	    for (auto xt : xState.in) {
	      const ProfileTransition& xTrans = x.trans[xt];
	      const XYCell& src = xyCell(xTrans.src,j);
	      log_accum_exp (imd, src(PairHMM::IMD) + xTrans.lpTrans);
	      log_accum_exp (iiw, src(PairHMM::IIW) + xTrans.lpTrans);
	    }
	}

	if (!yState.isNull()) {
	  // y-absorbing transitions into IDM, IMI
	  if (xState.isReady() || xEmpty) {
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
	  }
	  
	} else {  // yState.isNull()
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

	} else if (yState.isNull() && xState.isEmitOrStart()) {
	  // y-nonabsorbing transitions in IMM
	  for (auto yt : yState.in) {
	    const ProfileTransition& yTrans = y.trans[yt];
	    log_accum_exp (imm, cell(i,yTrans.src,PairHMM::IMM) + yTrans.lpTrans);
	  }

	} else {  // xState.isNull()
	  // x-nonabsorbing transitions in IMM
	  if (yState.isReady() || yEmpty)
	    for (auto xt : xState.in) {
	      const ProfileTransition& xTrans = x.trans[xt];
	      log_accum_exp (imm, cell(xTrans.src,j,PairHMM::IMM) + xTrans.lpTrans);
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
  auto clp = sourceTransitionsWithoutEmitOrAbsorb (destCell);

  const LogProb lpAbs = lpCellEmitOrAbsorb (destCell);
  for (auto& src_lp : clp)
    src_lp.second += lpAbs;

  return clp;
}

map<DPMatrix::CellCoords,LogProb> ForwardMatrix::sourceTransitionsWithoutEmitOrAbsorb (const CellCoords& destCell) {
  map<CellCoords,LogProb> clp;
  const ProfileState& xState = x.state[destCell.xpos];
  const ProfileState& yState = y.state[destCell.ypos];

  switch (destCell.state) {
  case PairHMM::IMD:
  case PairHMM::IIW:
    if (xState.isNull()) {
      // x-nonabsorbing transitions in IMD, IIW
      if (yState.isReady() || yEmpty)
	if (destCell.xpos < xSize - 1)
	  for (auto xt : xState.in)
	    clp[CellCoords(x.trans[xt].src,destCell.ypos,destCell.state)] = x.trans[xt].lpTrans;
    } else
      // x-absorbing transitions into IMD, IIW
      if (yState.isReady() || yEmpty)
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
      if (xState.isReady() || xEmpty)
	for (auto yt : yState.in)
	  for (auto s : hmm.sources (destCell.state))
	    clp[CellCoords(destCell.xpos,y.trans[yt].src,s)] = hmm.lpTrans(s,destCell.state) + y.trans[yt].lpTrans;
    break;

  case PairHMM::IMM:
    if (yState.isNull() && xState.isEmitOrStart()) {
      // y-nonabsorbing transitions in IMM
      if (destCell.ypos < ySize - 1)
	for (auto yt : yState.in)
	  clp[CellCoords(destCell.xpos,y.trans[yt].src,destCell.state)] = y.trans[yt].lpTrans;
    } else if (xState.isNull()) {
      // x-nonabsorbing transitions in IMM
      if (yState.isReady() || yEmpty)
	if (destCell.xpos < xSize - 1)
	  for (auto xt : xState.in)
	    clp[CellCoords(x.trans[xt].src,destCell.ypos,destCell.state)] = x.trans[xt].lpTrans;
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

string DPMatrix::toString (bool edgeOnly) const {
  ostringstream out;
  write (out, edgeOnly);
  return out.str();
}

void DPMatrix::write (ostream& out, bool edgeOnly) const {
  const auto states = PairHMM::states();
  CellCoords coords;
  for (coords.xpos = 0; coords.xpos < xSize - 1; ++coords.xpos)
    for (coords.ypos = 0; coords.ypos < ySize - 1; ++coords.ypos)
      if (edgeOnly ? atEdge(coords.xpos,coords.ypos) : inEnvelope(coords.xpos,coords.ypos))
	for (auto state: states) {
	  coords.state = state;
	  out << setw(16) << cell(coords)
	      << setw(6) << coords.xpos
	      << setw(6) << coords.ypos
	      << setw(6) << PairHMM::stateName(state,coords.xpos==0,coords.ypos==0)
	      << " " << x.tinyDescription(coords.xpos)
	      << " " << y.tinyDescription(coords.ypos)
	      << endl;
	}
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
  return (c.state == PairHMM::IMM && (x.state[c.xpos].isNull() || !y.state[c.ypos].isNull()))
    || c.state == PairHMM::IMD
    || c.state == PairHMM::IIW
    || c.state == PairHMM::EEE;
}

bool DPMatrix::changesY (const CellCoords& c) const {
  return (c.state == PairHMM::IMM && x.state[c.xpos].isEmitOrStart())
    || c.state == PairHMM::IDM
    || c.state == PairHMM::IMI
    || c.state == PairHMM::EEE;
}

list<DPMatrix::CellCoords> DPMatrix::equivAbsorbCells (const CellCoords& c) const {
  list<CellCoords> eq;
  if (c.state == PairHMM::IIW && !x.state[c.xpos].isNull())
    eq.push_back (CellCoords (c.xpos, c.ypos, PairHMM::IMD));
  else if (c.state == PairHMM::IMI && !y.state[c.ypos].isNull())
    eq.push_back (CellCoords (c.xpos, c.ypos, PairHMM::IDM));
  else if (changesX(c) && x.state[c.xpos].isNull() && x.equivAbsorbState.count(c.xpos))
    eq.push_back (CellCoords (x.equivAbsorbState.at(c.xpos), c.ypos, PairHMM::IMD));
  else if (changesY(c) && y.state[c.ypos].isNull() && y.equivAbsorbState.count(c.ypos))
    eq.push_back (CellCoords (c.xpos, y.equivAbsorbState.at(c.ypos), PairHMM::IDM));
  return eq;
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

ProfileState::SeqCoords ForwardMatrix::cellSeqCoords (const CellCoords& c) const {
  ProfileState::SeqCoords coords = x.state[c.xpos].seqCoords;
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
    else if (x.state[c.xpos].isEmitOrStart())
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
    if (!xNull && !yNull) {
      if (src.state == PairHMM::IMM || src.state == PairHMM::IMD) {
	c.indelCounts.insTime += hmm.l.t;
	c.indelCounts.delTime += hmm.l.t;
      }
      if (src.state == PairHMM::IMM || src.state == PairHMM::IDM) {
	c.indelCounts.insTime += hmm.r.t;
	c.indelCounts.delTime += hmm.r.t;
      }
    }
    break;
  case PairHMM::IMD:
    if (!xNull) {
      if (src.state == PairHMM::IMM || src.state == PairHMM::IMD) {
	c.indelCounts.insTime += hmm.l.t;
	c.indelCounts.delTime += hmm.l.t;
      }
      if (src.state == dest.state)
	c.indelCounts.delExt += 1;
      else {
	c.indelCounts.del += 1;
	c.indelCounts.delTime += hmm.r.delWait;
      }
    }
    break;
  case PairHMM::IIW:
    if (!xNull) {
      if (src.state == dest.state)
	c.indelCounts.insExt += 1;
      else {
	c.indelCounts.ins += 1;
	c.indelCounts.insTime += hmm.l.insWait;
      }
    }
    break;
  case PairHMM::IDM:
    if (!yNull) {
      if (src.state == PairHMM::IMM || src.state == PairHMM::IDM) {
	c.indelCounts.insTime += hmm.r.t;
	c.indelCounts.delTime += hmm.r.t;
      }
      if (src.state == dest.state)
	c.indelCounts.delExt += 1;
      else {
	c.indelCounts.del += 1;
	c.indelCounts.delTime += hmm.l.delWait;
      }
    }
    break;
  case PairHMM::IMI:
    if (!yNull) {
      if (src.state == dest.state)
	c.indelCounts.insExt += 1;
      else {
	c.indelCounts.ins += 1;
	c.indelCounts.insTime += hmm.r.insWait;
      }
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
  map<AlignRowIndex,SeqIdx> seqCoords;
  for (size_t n = 0; n < pv.size() - 1; ++n) {
    const AlignPath cap = cellAlignPath(pv[n]), tap = transitionAlignPath(pv[n],pv[n+1]);
    p = alignPathConcat (p, cap, tap);
    // do some consistency tests...
    for (const auto& rp: cap)
      seqCoords[rp.first] += alignPathResiduesInRow (rp.second);
    for (const auto& sc: x.state[pv[n].xpos].seqCoords)
      Assert (seqCoords[sc.first] == sc.second, "Sequence %d: cell x-coord is %d, path x-coord is %d", sc.first, sc.second, seqCoords[sc.first]);
    for (const auto& sc: y.state[pv[n].ypos].seqCoords)
      Assert (seqCoords[sc.first] == sc.second, "Sequence %d: cell y-coord is %d, path y-coord is %d", sc.first, sc.second, seqCoords[sc.first]);
    for (const auto& rp: tap)
      seqCoords[rp.first] += alignPathResiduesInRow (rp.second);
  }
  p = alignPathConcat (p, cellAlignPath(pv.back()));

  // check that parent & both children have entries (guard against bug that occurs if sequences are empty)
  ensureAlignPathHasRow (p, parentRowIndex);
  ensureAlignPathHasRow (p, x.rootRowIndex);
  ensureAlignPathHasRow (p, y.rootRowIndex);

  // report size
  const AlignRowIndex rows = p.size();
  const AlignColIndex cols = alignPathColumns (p);  // this will also test if alignment is flush
  LogThisAt(2,"Converted Forward matrix trace into alignment with " << rows << " rows and " << cols << " columns" << endl);

  return p;
}

Profile ForwardMatrix::makeProfile (const set<CellCoords>& cells, ProfilingStrategy strategy) {
  Profile prof (hmm.components(), alphSize, parentRowIndex);
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
	|| outgoingTransitionCount[c] > 1
	|| (strategy & KeepGapsOpen) != 0
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

  if (strategy & KeepGapsOpen)
    for (const auto& c : cells)
      if (!isAbsorbing(c) && profStateIndex.count(c)) {
	const auto equiv = equivAbsorbCells (c);
	if (equiv.size() && profStateIndex.count(equiv.front()))
	  prof.equivAbsorbState[profStateIndex[c]] = profStateIndex[equiv.front()];
      }

  if (strategy & CollapseChains)
    LogThisAt(5,"Eliminated " << plural(cells.size() - prof.state.size(), "state") << " from " << cells.size() << "-state profile using 'collapse chains' heuristic" << endl);
  
  // Calculate log-probabilities of "effective transitions" from cells to retained-cells.
  // Each effective transition represents a sum over paths.
  // A path is either a single direct transition (from source cell to destination retained-cell),
  // or a series of transitions starting from the source cell,
  // passing through one or more eliminated-cells, and stopping at the destination retained-cell.
  map<CellCoords,map<ProfileStateIndex,EffectiveTransition> > effTrans;  // effTrans[srcCell][destStateIdx]
  for (auto iter = cells.crbegin(); iter != cells.crend(); ++iter) {
    const CellCoords& iterCell = *iter;
    const map<CellCoords,LogProb>& slp = sourceTransitionsWithoutEmitOrAbsorb (iterCell);
    const LogProb cellLogProbInsert = eliminatedLogProbInsert (iterCell);
    if (profStateIndex.find(iterCell) != profStateIndex.end()) {
      // iterCell is to be retained. Incoming & outgoing paths can be kept separate
      const ProfileStateIndex cellIdx = profStateIndex[iterCell];
      for (const auto slpIter : slp) {
	const CellCoords& src = slpIter.first;
	const LogProb srcCellLogProbTrans = slpIter.second;
	EffectiveTransition& eff = effTrans[src][cellIdx];
	eff.lpPath = eff.lpBestAlignPath = srcCellLogProbTrans + cellLogProbInsert;
	eff.bestAlignPath = transitionAlignPath(src,iterCell);
	if (strategy & (CountSubstEvents | CountIndelEvents))
	  eff.counts = transitionEigenCounts(src,iterCell);
	// consistency check
	ProfileState::assertSeqCoordsConsistent (cellSeqCoords(src), prof.state[cellIdx], eff.bestAlignPath);
      }
    } else {
      // iterCell is to be eliminated. Connect incoming transitions & outgoing paths, summing iterCell out
      const auto& cellEffTrans = effTrans[iterCell];
      const AlignPath& cap = cellAlignPath (iterCell);
      EigenCounts cellCounts, srcCellCounts;
      if ((strategy & CountSubstEvents) != 0 && sumProd != NULL)
	cellCounts = cachedCellEigenCounts (iterCell, *sumProd);
      for (auto slpIter : slp) {
	const CellCoords& src = slpIter.first;
	const LogProb srcCellLogProbTrans = slpIter.second;
	if (strategy & (CountSubstEvents | CountIndelEvents))
	  srcCellCounts = transitionEigenCounts (src, iterCell) + cellCounts;
	auto& srcEffTrans = effTrans[src];
	for (const auto cellEffTransIter : cellEffTrans) {
	  const ProfileStateIndex& destIdx = cellEffTransIter.first;
	  const EffectiveTransition& cellDestEffTrans = cellEffTransIter.second;
	  EffectiveTransition& srcDestEffTrans = srcEffTrans[destIdx];
	  const LogProb lpPath = srcCellLogProbTrans + cellLogProbInsert + cellDestEffTrans.lpPath;
	  log_accum_exp (srcDestEffTrans.lpPath, lpPath);
	  if (strategy & (CountSubstEvents | CountIndelEvents)) {
	    const double ppPath = exp (lpPath - srcDestEffTrans.lpPath);
	    srcDestEffTrans.counts *= 1 - ppPath;
	    srcDestEffTrans.counts += (srcCellCounts + cellDestEffTrans.counts) * ppPath;
	  }
	  // we also want to keep the best single alignment consistent w/this set of paths
	  const LogProb srcDestLogProbBestAlignPath = srcCellLogProbTrans + cellLogProbInsert + cellDestEffTrans.lpBestAlignPath;
	  const AlignPath tap = transitionAlignPath(src,iterCell);
	  if (srcDestLogProbBestAlignPath > srcDestEffTrans.lpBestAlignPath) {
	    srcDestEffTrans.lpBestAlignPath = srcDestLogProbBestAlignPath;
	    srcDestEffTrans.bestAlignPath = alignPathConcat (tap, cap, cellDestEffTrans.bestAlignPath);
	  }
	  // consistency check
	  ProfileState::assertSeqCoordsConsistent (cellSeqCoords(iterCell), prof.state[destIdx], cellDestEffTrans.bestAlignPath);
	  ProfileState::assertSeqCoordsConsistent (cellSeqCoords(src), cellSeqCoords(iterCell), tap, cap);
	  ProfileState::assertSeqCoordsConsistent (cellSeqCoords(src), prof.state[destIdx], srcDestEffTrans.bestAlignPath);
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
      if (strategy & (CountSubstEvents | CountIndelEvents))
	trans.counts = srcDestEffTrans.counts;
      prof.trans.push_back (trans);
      (prof.state[destIdx].isNull() ? srcNullOut : srcAbsorbOut).push_back (transIdx);
      destIn.push_back (transIdx);
    }
  }

  prof.seq = x.seq;
  prof.seq.insert (y.seq.begin(), y.seq.end());

  // transform into ready/wait form & verify integrity
  prof.assertTransitionsConsistent();  // addReadyStates() will check this again
  prof.assertPathToEndExists();  // addReadyStates() will check this again
  prof = prof.addReadyStates();
  prof.assertSeqCoordsConsistent();
  
  return prof;
}

Profile ForwardMatrix::sampleProfile (random_engine& generator, size_t profileSamples, size_t maxCells, ProfilingStrategy strategy, size_t minLen, size_t maxLen) {
  map<CellCoords,size_t> cellCount;

  Require ((strategy & IncludeBestTrace) || profileSamples > 0, "Must allow at least one sample path in the profile");

  size_t nTraces = 0;
  if (strategy & IncludeBestTrace) {
    const Path best = bestTrace();
    for (auto& c : best)
      cellCount[c] = 2;  // avoid dropping these cells
    ++nTraces;
  }
  size_t nAccepted = 0;
  for (size_t n = 0; nAccepted < profileSamples && (maxCells == 0 || cellCount.size() < maxCells); ++n) {
    const Path sampled = sampleTrace (generator);
    if (LoggingThisAt(5)) {
      LogThisAt(5,"Trace #" << n+1 << ":");
      for (auto& c : sampled)
	LogThisAt(5," " << cellName(c));
      LogThisAt(5,endl);
    }
    size_t ancLen = 0;
    for (auto& c : sampled)
      switch (c.state) {
      case PairHMM::IMM:
      case PairHMM::IDM:
      case PairHMM::IMD:
	++ancLen;
      default:
	break;
      }
    if (ancLen < minLen || ancLen > maxLen)
      break;
    for (auto& c : sampled)
      ++cellCount[c];
    ++nTraces;
    ++nAccepted;
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
  EigenCounts c (hmm.components(), hmm.alphabetSize());
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
      } else if (x.state[cell.xpos].isEmitOrStart() && y.state[cell.ypos].isNull())
	col = y.alignColumn (cell.ypos);
      else if (x.state[cell.xpos].isNull())
	col = x.alignColumn (cell.xpos);
      break;
    case PairHMM::IMD:
      col = x.alignColumn (cell.xpos);
      if (!x.state[cell.xpos].isNull())
	col[parentRowIndex] = Alignment::wildcardChar;
      break;
    case PairHMM::IDM:
      col = y.alignColumn (cell.ypos);
      if (!y.state[cell.ypos].isNull())
	col[parentRowIndex] = Alignment::wildcardChar;
      break;
    case PairHMM::IIW:
      col = x.alignColumn (cell.xpos);
      break;
    case PairHMM::IMI:
      col = y.alignColumn (cell.ypos);
      break;
    default:
      break;
    }
  return col;
}

BackwardMatrix::BackwardMatrix (ForwardMatrix& fwd)
  : DPMatrix (fwd.x, fwd.y, fwd.hmm, fwd.envelope),
    fwd (fwd)
{
  lpEnd = 0;
  // transitions into EEE
  for (auto xt : x.end().in) {
    const ProfileTransition& xTrans = x.trans[xt];
    for (auto yt : y.end().in) {
      const ProfileTransition& yTrans = y.trans[yt];
      if (inEnvelope(xTrans.src,yTrans.src)) {
	XYCell& src = xyCell(xTrans.src,yTrans.src);
      
	src(PairHMM::IMM) = xTrans.lpTrans + yTrans.lpTrans + hmm.imm_eee;
	src(PairHMM::IMD) = xTrans.lpTrans + yTrans.lpTrans + hmm.imd_eee;
	src(PairHMM::IDM) = xTrans.lpTrans + yTrans.lpTrans + hmm.idm_eee;
	src(PairHMM::IMI) = xTrans.lpTrans + yTrans.lpTrans + hmm.imi_eee;
	src(PairHMM::IIW) = xTrans.lpTrans + yTrans.lpTrans + hmm.iiw_eee;
      }
    }
  }

  ProgressLog (plog, 5);
  plog.initProgress ("Backward algorithm (%s vs %s)", x.name.c_str(), y.name.c_str());

  for (int i = xSize - 2; i >= 0; --i) {
    const ProfileState& xState = x.state[i];

    plog.logProgress ((xSize - 2 - i) / (double) (xSize - 2), "state %d/%d", xSize - 1 - i, xSize);

    for (int j = ySize - 2; j >= 0; --j) {
      const ProfileState& yState = y.state[j];

      if (inEnvelope(i,j)) {

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
	if (yState.isReady() || yEmpty)
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
	if (xState.isReady() || xEmpty)
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
	if (yState.isReady() || yEmpty)
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
	  if (xState.isEmitOrStart())
	    log_accum_exp (imm, yTrans.lpTrans + dest(PairHMM::IMM));
	}
      }
    }
  }

  LogThisAt(6,"Backward log-likelihood is " << lpStart() << endl);
  if (gsl_fcmp (lpStart(), fwd.lpEnd, FWD_BACK_ERROR_TOLERANCE) != 0) {
    fwd.slowFillTest();
    slowFillTest();
    sourceDestTransTest();
    Warn ("Forward log-likelihood is %g, Backward log-likelihood is %g", fwd.lpEnd, lpStart());
  }
}

void ForwardMatrix::slowFillTest() {
  auto states = hmm.states();
  states.push_back (PairHMM::EEE);
  size_t nCells = 0, nTrans = 0;
  for (int i = 0; i < xSize; ++i)
    for (int j = 0; j < ySize; ++j)
      if (inEnvelope(i,j))
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
      if (inEnvelope(i,j))
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
      if (inEnvelope(i,j))
	for (auto s : states) {
	  const CellCoords testCell (i, j, s);
	  for (const auto& src_lp : fwd.sourceTransitions (testCell))
	    if (src_lp.second > -numeric_limits<double>::infinity()) {
	      auto destTrans = destTransitions (src_lp.first);
	      if (!destTrans.count(testCell))
		Warn ("Backward matrix is missing transition between %s and %s that is present in Forward matrix", cellName(src_lp.first).c_str(), cellName(testCell).c_str());
	      else
		Test (gsl_fcmp (src_lp.second, destTrans[testCell], FWD_BACK_ERROR_TOLERANCE) == 0, "Forward (%g) & Backward (%g) transitions between %s and %s don't match", src_lp.second, destTrans[testCell], cellName(src_lp.first).c_str(), cellName(testCell).c_str());
	    }
	  for (const auto& dest_lp : destTransitions (testCell))
	    if (dest_lp.second > -numeric_limits<double>::infinity()) {
	      auto srcTrans = fwd.sourceTransitions (dest_lp.first);
	      if (!srcTrans.count(testCell))
		Warn ("Forward matrix is missing transition between %s and %s that is present in Backward matrix", cellName(testCell).c_str(), cellName(dest_lp.first).c_str());
	      else
		Test (gsl_fcmp (dest_lp.second, srcTrans[testCell], FWD_BACK_ERROR_TOLERANCE) == 0, "Forward (%g) & Backward (%g) transitions between %s and %s don't match", srcTrans[testCell], dest_lp.second, cellName(testCell).c_str(), cellName(dest_lp.first).c_str());
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
  EigenCounts counts (hmm.components(), hmm.alphabetSize());
  counts.indelCounts.lp = fwd.lpEnd;
  
  const auto states = hmm.states();

  ProgressLog (plog, 4);
  plog.initProgress ("Forward-Backward counts (%s vs %s)", x.name.c_str(), y.name.c_str());

  for (ProfileStateIndex i = 0; i < xSize - 1; ++i) {
    const ProfileState& xState = x.state[i];
    plog.logProgress (i / (double) (xSize - 2), "state %d/%d", i + 1, xSize);

    for (ProfileStateIndex j = 0; j < ySize - 1; ++j) {
      const ProfileState& yState = y.state[j];

      if (inEnvelope(i,j)) {
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
  if (yState.isReady() || yEmpty)
    for (auto xt : xState.absorbOut) {
      const ProfileTransition& xTrans = x.trans[xt];
      clp[CellCoords(xTrans.dest,srcCell.ypos,PairHMM::IMD)] = hmm.lpTrans(srcCell.state,PairHMM::IMD) + xTrans.lpTrans;
      clp[CellCoords(xTrans.dest,srcCell.ypos,PairHMM::IIW)] = hmm.lpTrans(srcCell.state,PairHMM::IIW) + xTrans.lpTrans;
    }

  // y-absorbing transitions into IDM, IMI
  if (xState.isReady() || xEmpty)
    for (auto yt : yState.absorbOut) {
      const ProfileTransition& yTrans = y.trans[yt];
      clp[CellCoords(srcCell.xpos,yTrans.dest,PairHMM::IDM)] = hmm.lpTrans(srcCell.state,PairHMM::IDM) + yTrans.lpTrans;
      clp[CellCoords(srcCell.xpos,yTrans.dest,PairHMM::IMI)] = hmm.lpTrans(srcCell.state,PairHMM::IMI) + yTrans.lpTrans;
    }

  // x-nonabsorbing transitions in IMD, IIW, IMM
  if ((yState.isReady() || yEmpty) && (srcCell.state == PairHMM::IMD || srcCell.state == PairHMM::IIW || srcCell.state == PairHMM::IMM))
    for (auto xt : xState.nullOut) {
      const ProfileTransition& xTrans = x.trans[xt];
      if (xTrans.dest != xSize - 1)
	clp[CellCoords(xTrans.dest,srcCell.ypos,srcCell.state)] = xTrans.lpTrans;
    }

  // y-nonabsorbing transitions in IDM, IMI, IMM
  if (srcCell.state == PairHMM::IDM || srcCell.state == PairHMM::IMI || (xState.isEmitOrStart() && srcCell.state == PairHMM::IMM))
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

priority_queue<BackwardMatrix::CellPostProb> BackwardMatrix::cellsAbovePostProbThreshold (double minPostProb) const {
  priority_queue<CellPostProb> bc;
  const LogProb lppThreshold = log(minPostProb);
  const LogProb fwdEnd = fwd.lpEnd;
  const auto states = hmm.states();
  for (int i = xSize - 2; i >= 0; --i)
    for (int j = ySize - 2; j >= 0; --j)
      if (inEnvelope(i,j)) {
	const XYCell& backSrc = xyCell(i,j);
	const XYCell& fwdSrc = fwd.xyCell(i,j);
	for (auto s : states) {
	  const LogProb lpp = backSrc(s) + fwdSrc(s) - fwdEnd;
	  if (lpp >= lppThreshold)
	    bc.push (CellPostProb (i, j, s, lpp));
	}
      }
  return bc;
}

Profile BackwardMatrix::bestProfile (ProfilingStrategy strategy) {
  set<CellCoords> cells;
  addTrace (endCell, cells, 0, (strategy & KeepGapsOpen) != 0);
  return fwd.makeProfile (cells, strategy);
}

Profile BackwardMatrix::postProbProfile (double minPostProb, size_t maxCells, ProfilingStrategy strategy) {
  priority_queue<CellPostProb> bc = cellsAbovePostProbThreshold (minPostProb);
  set<CellCoords> cells;
  if (bc.empty() || (strategy & IncludeBestTrace))
    addCells (cells, 0, fwd.bestTrace(), list<CellCoords>(), (strategy & KeepGapsOpen) != 0);
  while ((maxCells == 0 || cells.size() < maxCells) && !bc.empty()) {
    const CellCoords& best = bc.top();
    if (cells.count (best))
      bc.pop();
    else
      if (!addTrace (best, cells, maxCells, (strategy & KeepGapsOpen) != 0))
	break;
  }
  return fwd.makeProfile (cells, strategy);
}

bool BackwardMatrix::addCells (set<CellCoords>& cells, size_t maxCells, const list<CellCoords>& fwdTrace, const list<CellCoords>& backTrace, bool keepGapsOpen) {
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
    return false;
  if (LoggingThisAt(6)) {
    LogThisAt(6,"Adding the following cells to profile:");
    for (const auto& c : newCells)
      LogThisAt(6," " << cellName(c));
    LogThisAt(6,endl);
  }
  cells.insert (newCells.begin(), newCells.end());
  if (keepGapsOpen)
    for (const auto& newCell : newCells) {
      const list<CellCoords>& eqvCells = equivAbsorbCells (newCell);
      for (auto& eqvCell : eqvCells)
	if (!cells.count (eqvCell) && cellPostProb(eqvCell) > 0 && inEnvelope(eqvCell.xpos,eqvCell.ypos))
	  addTrace (eqvCell, cells, maxCells, false);
    }
  return true;
}

bool BackwardMatrix::addTrace (const CellCoords& cell, set<CellCoords>& cells, size_t maxCells, bool keepGapsOpen) {
  LogThisAt(5,"Starting traceback/forward from " << cellName(cell) << endl);
  const list<CellCoords> fwdTrace = fwd.bestTrace(cell), backTrace = bestTrace(cell);
  return addCells (cells, maxCells, fwdTrace, backTrace, keepGapsOpen);
  return true;
}
