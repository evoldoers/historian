#ifndef QUICKALIGN_INCLUDED
#define QUICKALIGN_INCLUDED

#include "diagenv.h"
#include "model.h"
#include "alignpath.h"

class QuickAlignMatrix {
public:
  enum State { Start, Match, Insert, Delete };
  const DiagonalEnvelope* penv;
  const FastSeq *px, *py;
  vguard<AlphTok> xTok, yTok;
  SeqIdx xLen, yLen;
  vguard<LogProb> cell;
  LogProb start, end, result;
  static double dummy;

  const RateModel& model;
  const double time;
  vguard<vguard<LogProb> > submat;  // log odds-ratio
  LogProb m2m, m2i, m2d, m2e, i2i, i2m, i2d, i2e, d2d, d2m, d2e;
  
  QuickAlignMatrix (const DiagonalEnvelope& env, const RateModel& model, double time);
  inline LogProb& getCell (SeqIdx i, SeqIdx j, unsigned int offset) {
    const int storageIndex = penv->getStorageIndexUnsafe (i, j);
    return cell[storageIndex*3 + offset];
  }
  inline const LogProb& getCell (SeqIdx i, SeqIdx j, unsigned int offset) const {
    const int storageIndex = penv->getStorageIndexSafe (i, j);
    return storageIndex < 0 ? dummy : cell[storageIndex*3 + offset];
  }
  inline LogProb& mat (SeqIdx i, SeqIdx j) { return getCell(i,j,0); }
  inline LogProb& ins (SeqIdx i, SeqIdx j) { return getCell(i,j,1); }
  inline LogProb& del (SeqIdx i, SeqIdx j) { return getCell(i,j,2); }
  inline const LogProb& mat (SeqIdx i, SeqIdx j) const { return getCell(i,j,0); }
  inline const LogProb& ins (SeqIdx i, SeqIdx j) const { return getCell(i,j,1); }
  inline const LogProb& del (SeqIdx i, SeqIdx j) const { return getCell(i,j,2); }

  inline double matchEmitScore (SeqIdx i, SeqIdx j) const {
    Assert (i > 0 && j > 0 && i <= xLen && j <= yLen, "Out of range: (i,j)=(%u,%u) (xLen,yLen)=(%u,%u)", i, j, xLen, yLen);
    return submat[xTok[i-1]][yTok[j-1]];
  }

  LogProb cellScore (SeqIdx i, SeqIdx j, State state) const;
  static const char* stateToString (State state);
  static size_t cellSize() { return 3*sizeof(double); }

  bool resultIsFinite() const { return result > -numeric_limits<double>::infinity(); }
  AlignPath alignment() const;

protected:
  static void updateMax (LogProb& currentMax, State& currentMaxIdx, double candidateMax, State candidateMaxIdx);
};



#endif /* QUICKALIGN_INCLUDED */
