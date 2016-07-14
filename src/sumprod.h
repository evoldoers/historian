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

struct SumProductStorage {
  // F_n(x_n): variable->function, tip->root messages
  // E_n(x_p): function->variable, tip->root messages
  // G_n(x_n): function->variable, root->tip messages
  // G_p(x_p)*E_s(x_p): variable->function, root->tip messages
  vguard<vguard<vguard<double> > > E, F, G;
  vguard<vguard<LogProb> > logE, logF, logG;  // logs of rescaling factors, used to prevent underflow
  
  vguard<char> gappedCol;
  vguard<AlignRowIndex> ungappedRows, roots;

  vguard<LogProb> cptLogLike;
  LogProb colLogLike;  // marginal likelihood, all unobserved states summed out

  SumProductStorage (size_t components, size_t nodes, size_t alphabetSize);
  SumProductStorage() { }
};

class SumProduct : private SumProductStorage {
private:
  void assertSingleRoot() const;
  
public:
  const RateModel& model;
  const Tree& tree;

  vguard<TreeNodeIndex> preorder, postorder;  // modify these to visit only subsets of nodes. postorder for fillUp(), preorder for fillDown()
  
  vguard<vguard<double> > insProb;  // insProb[cpt][state]
  vguard<vguard<vguard<vguard<double> > > > branchSubProb;  // branchSubProb[cpt][node][parentState][nodeState]

  EigenModel eigen;
  vguard<gsl_matrix_complex*> branchEigenSubCount;
  
  SumProduct (const RateModel& model, const Tree& tree);
  ~SumProduct();

  void initColumn (const map<AlignRowIndex,char>& seq);
  AlignRowIndex columnRoot() const;

  inline const vguard<AlignRowIndex>& ungappedRowIndices() const { return ungappedRows; }
  inline LogProb columnLogLikelihood() const { return colLogLike; }
  
  inline bool isGap (AlignRowIndex row) const { return Alignment::isGap (gappedCol[row]); }
  inline bool isWild (AlignRowIndex row) const { return Alignment::isWildcard (gappedCol[row]); }
  inline bool columnEmpty() const { return ungappedRows.empty(); }

  LogProb computeColumnLogLikelihoodAt (AlignRowIndex row) const;

  void fillUp();  // E, F
  void fillDown();  // G
  
  vguard<LogProb> logNodePostProb (AlignRowIndex node) const;  // marginalizes out component
  vguard<vguard<LogProb> > logNodeExcludedPostProb (TreeNodeIndex node, TreeNodeIndex exclude, bool normalize = true) const;
  LogProb logBranchPostProb (int component, AlignRowIndex node, AlphTok parentState, AlphTok nodeState) const;
  AlphTok maxPostState (AlignRowIndex node) const;  // maximum a posteriori reconstruction, component marginalized

  void accumulateEigenCounts (vguard<vguard<double> >& rootCounts, vguard<vguard<vguard<gsl_complex> > >& eigenCounts, double weight = 1.) const;
  void accumulateSubCounts (vguard<vguard<double> >& rootCounts, vguard<vguard<vguard<double> > >& subCounts, double weight = 1) const;

private:
  void initColumn();  // populates ungappedRows
  void accumulateRootCounts (vguard<vguard<double> >& rootCounts, double weight = 1) const;
  
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
