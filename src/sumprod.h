#ifndef SUMPROD_INCLUDED
#define SUMPROD_INCLUDED

#include <gsl/gsl_matrix_complex_double.h>
#include <gsl/gsl_vector_complex_double.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

#include "model.h"
#include "tree.h"
#include "fastseq.h"
#include "alignpath.h"

class EigenModel {
public:
  const RateModel& model;

  gsl_vector_complex *eval;
  gsl_matrix_complex *evec;  // right eigenvectors
  gsl_matrix_complex *evecInv;  // left eigenvectors

  EigenModel (const RateModel& model);
  ~EigenModel();

  gsl_matrix_complex* eigenSubCount (double t) const;
  double getSubCount (AlphTok a, AlphTok b, AlphTok i, AlphTok j, const gsl_matrix* sub, const gsl_matrix_complex* eSubCount) const;
  void accumSubCounts (vguard<vguard<double> >& count, AlphTok a, AlphTok b, double weight, const gsl_matrix* sub, const gsl_matrix_complex* eSubCount) const;

  // for testing purposes...
  double getSubProb (double t, AlphTok i, AlphTok j) const;
  gsl_matrix* getSubProbMatrix (double t) const;
  gsl_matrix_complex* getRateMatrix() const;
  gsl_matrix_complex* evecInv_evec() const;

  vguard<vguard<double> > getSubCounts (const vguard<vguard<gsl_complex> >& eigenCounts) const;
  
private:
  vguard<gsl_complex> ev, ev_t, exp_ev_t;
  void compute_exp_ev_t (double t);
  double getSubProbInner (double t, AlphTok i, AlphTok j) const;
  
  EigenModel (const EigenModel&) = delete;
  EigenModel& operator= (const EigenModel&) = delete;
};

class SumProduct {
public:
  const RateModel& model;
  const Tree& tree;

  vguard<TreeNodeIndex> preorder, postorder;
  
  vguard<LogProb> logInsProb;
  vguard<vguard<vguard<LogProb> > > branchLogSubProb;  // branchLogSubProb[node][parentState][nodeState]

  EigenModel eigen;
  vguard<gsl_matrix_complex*> branchEigenSubCount;

  vguard<char> gappedCol;
  vguard<AlignRowIndex> ungappedRows;

  // F_n(x_n): variable->function, tip->root messages
  // E_n(x_p): function->variable, tip->root messages
  // G_n(x_n): function->variable, root->tip messages
  // G_p(x_p)*E_s(x_p): variable->function, root->tip messages
  vguard<vguard<LogProb> > logE, logF, logG;
  LogProb colLogLike;  // marginal likelihood, all unobserved states summed out
  
  SumProduct (const RateModel& model, const Tree& tree);
  ~SumProduct();

  void initColumn (const map<AlignRowIndex,char>& seq);
  
  inline bool isGap (AlignRowIndex row) const { return Alignment::isGap (gappedCol[row]); }
  inline bool isWild (AlignRowIndex row) const { return Alignment::isWildcard (gappedCol[row]); }
  inline bool columnEmpty() const { return ungappedRows.empty(); }
  inline AlignRowIndex root() const { return ungappedRows.back(); }

  void fillUp();  // E, F
  void fillDown();  // G

  vguard<LogProb> logNodePostProb (AlignRowIndex node) const;
  vguard<LogProb> logNodeExcludedPostProb (AlignRowIndex node, AlignRowIndex exclude) const;
  LogProb logBranchPostProb (AlignRowIndex node, AlphTok parentState, AlphTok nodeState) const;
  AlphTok maxPostState (AlignRowIndex node) const;  // maximum a posteriori reconstruction

  void accumulateEigenCounts (vguard<double>& rootCounts, vguard<vguard<gsl_complex> >& eigenCounts, double weight = 1.) const;
  vguard<vguard<double> > getSubCounts (vguard<vguard<gsl_complex> >& eigenCounts) const;  // wait times on diagonal

  void accumulateSubCounts (vguard<double>& rootCounts, vguard<vguard<double> >& subCounts, double weight = 1) const;

private:
  void initColumn();  // populates ungappedRows
  void accumulateRootCounts (vguard<double>& rootCounts, double weight = 1) const;
  
  SumProduct (const SumProduct&) = delete;
  SumProduct& operator= (const SumProduct&) = delete;
};

class AlignColSumProduct : public SumProduct {
public:
  typedef map<AlignRowIndex,map<AlignColIndex,map<char,double> > > ReconPostProbMap;

  const vguard<FastSeq>& gapped;  // tree node index must match alignment row index
  AlignColIndex col;
  
  AlignColSumProduct (const RateModel& model, const Tree& tree, const vguard<FastSeq>& gapped);

  bool alignmentDone() const;
  void nextColumn();

  void appendAncestralReconstructedColumn (vguard<FastSeq>& out) const;
  void appendAncestralPostProbColumn (ReconPostProbMap& out, double minProb = .01, double maxProb = .999) const;
  
private:
  void initAlignColumn();  // populates ungappedRows
};

#endif /* SUMPROD_INCLUDED */
