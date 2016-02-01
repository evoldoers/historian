#ifndef SUMPROD_INCLUDED
#define SUMPROD_INCLUDED

#include <gsl/gsl_matrix_complex_double.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

#include "model.h"
#include "tree.h"
#include "fastseq.h"
#include "alignpath.h"

class AlignColSumProduct {
public:
  const RateModel& model;
  const Tree& tree;
  const vguard<FastSeq>& gapped;

  vguard<LogProb> logInsProb;
  vguard<vguard<vguard<LogProb> > > branchLogSubProb;  // branchLogSubProb[node][src][dest]

  gsl_vector_complex *eval;
  gsl_matrix_complex *evec;  // right eigenvectors
  gsl_matrix_complex *evecTrans;  // left eigenvectors

  AlignColIndex col;
  vguard<SeqIdx> seqPos;
  vguard<AlignRowIndex> ungappedRows;

  // F_n(x_n): variable->function, tip->root messages
  // E_n(x_p): function->variable, tip->root messages
  // G_n(x_n): function->variable, root->tip messages
  // G_p(x_p)*E_s(x_p): variable->function, root->tip messages
  vguard<LogProb> logE, logF, logG;

  AlignColSumProduct (const RateModel& model, const Tree& tree, const vguard<FastSeq>& gapped);
  ~AlignColSumProduct();

  bool alignmentDone() const;
  void nextColumn();

  inline bool columnEmpty() const { return ungappedRows.empty(); }
  inline AlignRowIndex root() const { return ungappedRows.back(); }

  void fillUp();  // E, F
  void fillDown();  // G

  LogProb logColumnProb() const;
  vguard<LogProb> logNodePostProb (AlignRowIndex node) const;
  AlphTok maxPostState (AlignRowIndex node) const;  // maximum a posteriori reconstruction

  void accumulateCounts (gsl_vector* rootCounts, gsl_matrix_complex* eigenCounts) const;
  gsl_matrix* getSubCounts (gsl_matrix_complex* eigenCounts) const;

private:
  void initColumn();  // populates ungappedRows

  AlignColSumProduct (const AlignColSumProduct&) = delete;
  AlignColSumProduct& operator= (const AlignColSumProduct&) = delete;
};

#endif /* SUMPROD_INCLUDED */
