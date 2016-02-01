#include <gsl/gsl_eigen.h>
#include "sumprod.h"
#include "util.h"

AlignColSumProduct::AlignColSumProduct (const RateModel& model, const Tree& tree, const vguard<FastSeq>& gapped)
  : model (model),
    tree (tree),
    gapped (gapped),
    logInsProb (log_gsl_vector (model.insProb)),
    branchLogSubProb (tree.nodes()),
    eval (gsl_vector_complex_alloc (model.alphabetSize())),
    evec (gsl_matrix_complex_alloc (model.alphabetSize(), model.alphabetSize())),
    col (0),
    logE (tree.nodes(), vguard<LogProb> (model.alphabetSize())),
    logF (tree.nodes(), vguard<LogProb> (model.alphabetSize())),
    logG (tree.nodes(), vguard<LogProb> (model.alphabetSize()))
{
  Require (tree.nodes() == gapped.size(), "Every tree node must have an alignment row");

  for (AlignRowIndex r = 0; r < gapped.size() - 1; ++r) {
    ProbModel pm (model, tree.branchLength(r));
    LogProbModel lpm (pm);
    branchLogSubProb[r].swap (lpm.logSubProb);
  }

  initColumn();

  gsl_eigen_nonsymmv_workspace *workspace = gsl_eigen_nonsymmv_alloc (model.alphabetSize());
  CheckGsl (gsl_eigen_nonsymmv (model.subRate, eval, evec, workspace));
  gsl_eigen_nonsymmv_free (workspace);
}

AlignColSumProduct::~AlignColSumProduct() {
  gsl_vector_complex_free (eval);
  gsl_matrix_complex_free (evec);
}

void AlignColSumProduct::initColumn() {
  ungappedRows.clear();
  vguard<int> ungappedKids (tree.nodes(), 0);
  vguard<TreeNodeIndex> roots;
  for (TreeNodeIndex r = 0; r < tree.nodes(); ++r)
    if (!isGap(r)) {
      ungappedRows.push_back (r);
      Assert (isWild(r) || ungappedKids[r] == 0, "Internal node sequences must be wildcards (%c)", Alignment::gapChar);
      const TreeNodeIndex rp = tree.parentNode(r);
      if (rp < 0 || isGap(rp))
	roots.push_back (r);
      else
	++ungappedKids[rp];
    }
  Assert (roots.size() == 1, "Multiple root nodes");
}

bool AlignColSumProduct::alignmentDone() const {
  return col >= gapped.front().length();
}

void AlignColSumProduct::nextColumn() {
  ++col;
  initColumn();
}

void AlignColSumProduct::fillUp() {
  for (auto r : ungappedRows) {
    if (isWild(r))
      for (AlphTok i = 0; i < model.alphabetSize(); ++i)
	logF[i] = 0;
    else {
      for (AlphTok i = 0; i < model.alphabetSize(); ++i)
	logF[i] = -numeric_limits<double>::infinity();
      logF[model.tokenize(gapped[r].seq[col])] = 0;
    }

    if (r == root())
      logLike = logInnerProduct (logF[r], logInsProb);
    else {
      // WRITE ME
    }
  }
}
