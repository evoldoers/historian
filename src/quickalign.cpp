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
    xTok (env.px->tokens (model.alphabet)),
    yTok (env.py->tokens (model.alphabet)),
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
  submat = lpm.logSubProb;
  for (AlphTok i = 0; i < model.alphabetSize(); ++i)
    for (AlphTok j = 0; j < model.alphabetSize(); ++j)
      submat[i][j] -= lpm.logInsProb[j];

  const double gapProb = pm.ins + (1 - pm.ins) * pm.del;
  const double noGapProb = 1 - gapProb;
  const double gapExt = 1 / ( (pm.ins/gapProb) / pm.insExt + (1 - pm.ins/gapProb) / pm.delExt );
  const double noGapExt = 1 - gapExt;

  m2i = log (gapProb);
  m2d = log (noGapProb * gapProb);
  m2m = log (noGapProb * noGapProb);
  m2e = 0;

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

  ProgressLog (plog, 4);
  plog.initProgress ("Viterbi algorithm (%s vs %s)", x.name.c_str(), y.name.c_str());

  start = 0;
  for (SeqIdx j = 1; j <= yLen; ++j) {

    plog.logProgress (j / (double) yLen, "base %d/%d", j, yLen);

    for (DiagonalEnvelope::iterator pi = env.begin(j);
	 !pi.finished();
	 ++pi) {

      const SeqIdx i = *pi;

      mat(i,j) = max (max (mat(i-1,j-1) + m2m,
			   del(i-1,j-1) + d2m),
		      ins(i-1,j-1) + i2m);

      if (j == 1 && i == 1)
	mat(i,j) = max (mat(i,j),
			start);

      mat(i,j) += matchEmitScore(i,j);

      ins(i,j) = max (ins(i,j-1) + i2i,
		      mat(i,j-1) + m2i);

      del(i,j) = max (max (ins(i-1,j) + i2d,
			   del(i-1,j) + d2d),
		      mat(i-1,j) + m2d);

      if (j == yLen && i == xLen)
	end = max (max (max (end,
			     mat(i,j) + m2e),
			ins(i,j) + i2e),
		   del(i,j) + d2e);
    }
  }

  result = end;

  LogThisAt(4, "Viterbi score: " << result << endl);
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

AlignPath QuickAlignMatrix::alignment() const {
  Require (resultIsFinite(), "Can't do Viterbi traceback if final score is -infinity");
  AlignPath path;
  SeqIdx i = xLen, j = yLen;
  State state = Match;
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
      if (j == 0 && i == 0)
	updateMax (srcSc, state, emitSc, Start);
      Assert (srcSc == mat(i+1,j+1), "Traceback error");
      break;

    case Insert:
      --j;
      path[0].insert (path[0].begin(), false);
      path[1].insert (path[1].begin(), true);
      updateMax (srcSc, state, mat(i,j) + m2i, Match);
      updateMax (srcSc, state, ins(i,j) + i2i, Insert);
      Assert (srcSc == ins(i,j+1), "Traceback error");
      break;

    case Delete:
      --i;
      path[0].insert (path[0].begin(), true);
      path[1].insert (path[1].begin(), false);
      updateMax (srcSc, state, mat(i,j) + m2d, Match);
      updateMax (srcSc, state, ins(i,j) + i2d, Insert);
      updateMax (srcSc, state, del(i,j) + d2d, Delete);
      Assert (srcSc == del(i+1,j), "Traceback error");
      break;

    default:
      Abort ("Traceback error");
      break;
    }
  }
  return path;
}
