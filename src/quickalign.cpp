#include <math.h>
#include "quickalign.h"
#include "logger.h"

LogProb QuickAlignMatrix::dummy = -numeric_limits<double>::infinity();

QuickAlignMatrix::QuickAlignMatrix (const DiagonalEnvelope& env, const RateModel& model, double time)
  : penv (&env),
    px (env.px),
    py (env.py),
    xLen (px->length()),
    yLen (py->length()),
    xTok (env.px->unvalidatedTokens (model.alphabet)),
    yTok (env.py->unvalidatedTokens (model.alphabet)),
    cell (env.totalStorageSize * 3, -numeric_limits<double>::infinity()),
    start (-numeric_limits<double>::infinity()),
    end (-numeric_limits<double>::infinity()),
    result (-numeric_limits<double>::infinity()),
    model (model),
    time (time)
{
  // compute scores
  ProbModel pm (model, time);
  LogProbModel lpm (pm);
  submat = lpm.logSubProb.front();
  for (AlphTok i = 0; i < model.alphabetSize(); ++i)
    for (AlphTok j = 0; j < model.alphabetSize(); ++j)
      submat[i][j] -= lpm.logInsProb.front()[j];

  const double gapProb = pm.ins + (1 - pm.ins) * pm.del;
  const double noGapProb = 1 - gapProb;
  const double gapExt = 1 / ( (pm.ins/gapProb) / pm.insExt + (1 - pm.ins/gapProb) / pm.delExt );
  const double noGapExt = 1 - gapExt;

  noGap = log(noGapProb);
  gapOpen = log(gapProb) + log(noGapExt);
  gapExtend = log(gapExt);
  
  m2i = log (gapProb);
  m2d = log (noGapProb * gapProb);
  m2m = log (noGapProb * noGapProb);

  i2i = log (gapExt);
  i2d = log (noGapExt * gapProb);
  i2m = log (noGapExt * noGapProb);
  i2e = i2m;

  d2d = log (gapExt);
  d2m = log (noGapExt);
  d2e = d2m;
  
  // fill
  const FastSeq& x (*px);
  const FastSeq& y (*py);

  ProgressLog (plog, 5);
  plog.initProgress ("Viterbi algorithm (%s vs %s)", x.name.c_str(), y.name.c_str());

  start = 0;
  xEnd = yEnd = 0;
  for (SeqIdx j = 1; j <= yLen; ++j) {

    plog.logProgress (j / (double) yLen, "base %d/%d", j, yLen);

    for (DiagonalEnvelope::iterator pi = env.begin(j);
	 !pi.finished();
	 ++pi) {

      const SeqIdx i = *pi;

      mat(i,j) = max (max (mat(i-1,j-1) + m2m,
			   del(i-1,j-1) + d2m),
		      ins(i-1,j-1) + i2m);

      mat(i,j) = max (mat(i,j),
		      start + startGapScore(i,j));

      mat(i,j) += matchEmitScore(i,j);

      ins(i,j) = max (ins(i,j-1) + i2i,
		      mat(i,j-1) + m2i);

      del(i,j) = max (max (ins(i-1,j) + i2d,
			   del(i-1,j) + d2d),
		      mat(i-1,j) + m2d);

      const LogProb ijEnd = mat(i,j) + endGapScore(i,j);
      if (ijEnd > end) {
	xEnd = i;
	yEnd = j;
	end = ijEnd;
      }
    }
  }

  result = end;

  LogThisAt(6, "Viterbi score: " << result << endl);
}

LogProb QuickAlignMatrix::cellScore (SeqIdx i, SeqIdx j, State state) const {
  LogProb cs = numeric_limits<double>::quiet_NaN();
  switch (state) {
  case Match:
    cs = mat(i,j);
    break;
  case Insert:
    cs = ins(i,j);
    break;
  case Delete:
    cs = del(i,j);
    break;
  default:
    break;
  }
  return cs;
}

const char* QuickAlignMatrix::stateToString (State state) {
  const char* s = "Unknown";
  switch (state) {
  case Start:
    s = "Start";
    break;
  case Match:
    s = "Match";
    break;
  case Insert:
    s = "Insert";
    break;
  case Delete:
    s = "Delete";
    break;
  default:
    break;
  }
  return s;
}

void QuickAlignMatrix::updateMax (double& currentMax, State& currentMaxIdx, double candidateMax, State candidateMaxIdx) {
  if (candidateMax > currentMax) {
    currentMax = candidateMax;
    currentMaxIdx = candidateMaxIdx;
  }
}

AlignPath QuickAlignMatrix::alignPath() const {
  Require (resultIsFinite(), "Can't do Viterbi traceback if final score is -infinity");
  SeqIdx i = xEnd, j = yEnd;
  State state = Match;
  Assert (i > 0 && j > 0, "Traceback error at (%lu,%lu,End)", i, j);
  AlignPath path;
  path[0] = vector<bool> (xLen - xEnd, true);
  path[1] = vector<bool> (xLen - xEnd, false);
  path[0].insert (path[0].end(), yLen - yEnd, false);
  path[1].insert (path[1].end(), yLen - yEnd, true);
  while (state != Start) {
    LogThisAt(9, "Traceback: i=" << i << " j=" << j << " state=" << stateToString(state) << " score=" << cellScore(i,j,state) << endl);
    LogProb srcSc = -numeric_limits<double>::infinity();
    LogProb emitSc = 0;
    switch (state) {
    case Match:
      emitSc = matchEmitScore(i,j);
      --i;
      --j;
      path[0].insert (path[0].begin(), true);
      path[1].insert (path[1].begin(), true);
      updateMax (srcSc, state, mat(i,j) + m2m + emitSc, Match);
      updateMax (srcSc, state, ins(i,j) + i2m + emitSc, Insert);
      updateMax (srcSc, state, del(i,j) + d2m + emitSc, Delete);
      updateMax (srcSc, state, start + startGapScore(i+1,j+1) + emitSc, Start);
      Assert (srcSc == mat(i+1,j+1), "Traceback error at (%lu,%lu,Match)", i+1, j+1);
      break;

    case Insert:
      --j;
      path[0].insert (path[0].begin(), false);
      path[1].insert (path[1].begin(), true);
      updateMax (srcSc, state, mat(i,j) + m2i, Match);
      updateMax (srcSc, state, ins(i,j) + i2i, Insert);
      Assert (srcSc == ins(i,j+1), "Traceback error at (%lu,%lu,Insert)", i, j+1);
      break;

    case Delete:
      --i;
      path[0].insert (path[0].begin(), true);
      path[1].insert (path[1].begin(), false);
      updateMax (srcSc, state, mat(i,j) + m2d, Match);
      updateMax (srcSc, state, ins(i,j) + i2d, Insert);
      updateMax (srcSc, state, del(i,j) + d2d, Delete);
      Assert (srcSc == del(i+1,j), "Traceback error at (%lu,%lu,Delete)", i+1, j);
      break;

    default:
      Abort ("Traceback error");
      break;
    }
  }
  path[0].insert (path[0].begin(), i, true);
  path[1].insert (path[1].begin(), i, false);
  path[0].insert (path[0].begin(), j, false);
  path[1].insert (path[1].begin(), j, true);
  LogThisAt(8,"Traceback alignment has " << plural(alignPathColumns(path),"column") << endl);  // checks alignment is flush
  Assert (alignPathResiduesInRow (path[0]) == xLen, "Traceback error: x row has %lu steps, expected %lu", alignPathResiduesInRow (path[0]), xLen);
  Assert (alignPathResiduesInRow (path[1]) == yLen, "Traceback error: y row has %lu steps, expected %lu", alignPathResiduesInRow (path[1]), yLen);
  return path;
}

AlignPath QuickAlignMatrix::alignPath (AlignRowIndex row1, AlignRowIndex row2) const {
  AlignPath oldPath = alignPath();
  AlignPath newPath;
  newPath[row1] = oldPath[0];
  newPath[row2] = oldPath[1];
  return newPath;
}

Alignment QuickAlignMatrix::alignment() const {
  AlignPath path = alignPath();
  vguard<FastSeq> seqs;
  seqs.push_back (*px);
  seqs.push_back (*py);
  return Alignment (seqs, path);
}

vguard<FastSeq> QuickAlignMatrix::gappedSeq() const {
  return alignment().gapped();
}
