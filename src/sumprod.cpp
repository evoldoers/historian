#include <gsl/gsl_eigen.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_math.h>

#include "sumprod.h"
#include "util.h"
#include "logger.h"

#define SUMPROD_RESCALE_THRESHOLD 1e-30

#define SUMPROD_EPSILON 1e-6
#define SUMPROD_NEAR_EQ(X,Y) (gsl_fcmp (X, Y, SUMPROD_EPSILON) == 0)
#define SUMPROD_NEAR_EQ_COMPLEX(X,Y) (SUMPROD_NEAR_EQ(GSL_REAL(X),GSL_REAL(Y)) && SUMPROD_NEAR_EQ(GSL_IMAG(X),GSL_IMAG(Y)))
#define SUMPROD_NEAR_REAL(X) (abs(GSL_IMAG(X)) < SUMPROD_EPSILON)

string tempComplexMatrixToString (gsl_matrix_complex* mx);
string complexMatrixToString (const gsl_matrix_complex* mx);
string complexVectorToString (const gsl_vector_complex* v);

string complexMatrixToString (const vguard<vguard<gsl_complex> >& mx);
string complexVectorToString (const vector<gsl_complex>& v);

EigenModel::EigenModel (const RateModel& model)
  : model (model),
    eval (gsl_vector_complex_alloc (model.alphabetSize())),
    evec (gsl_matrix_complex_alloc (model.alphabetSize(), model.alphabetSize())),
    evecInv (gsl_matrix_complex_alloc (model.alphabetSize(), model.alphabetSize())),
    ev (model.alphabetSize()),
    ev_t (model.alphabetSize()),
    exp_ev_t (model.alphabetSize())
{
  gsl_matrix *R = gsl_matrix_alloc (model.alphabetSize(), model.alphabetSize());
  gsl_matrix_memcpy (R, model.subRate);
  
  gsl_eigen_nonsymmv_workspace *workspace = gsl_eigen_nonsymmv_alloc (model.alphabetSize());
  CheckGsl (gsl_eigen_nonsymmv (R, eval, evec, workspace));
  gsl_eigen_nonsymmv_free (workspace);
  gsl_matrix_free (R);

  gsl_matrix_complex *LU = gsl_matrix_complex_alloc (model.alphabetSize(), model.alphabetSize());
  gsl_permutation *perm = gsl_permutation_alloc (model.alphabetSize());
  int permSig = 0;
  gsl_matrix_complex_memcpy (LU, evec);
  CheckGsl (gsl_linalg_complex_LU_decomp (LU, perm, &permSig));
  CheckGsl (gsl_linalg_complex_LU_invert (LU, perm, evecInv));
  gsl_matrix_complex_free (LU);
  gsl_permutation_free (perm);

  for (AlphTok i = 0; i < model.alphabetSize(); ++i)
    ev[i] = gsl_vector_complex_get (eval, i);

  LogThisAt(8,"Eigenvalues:" << complexVectorToString(ev) << endl
	    << "Right eigenvector matrix, V:" << endl << complexMatrixToString(evec)
	    << "Left eigenvector matrix, V^{-1}:" << endl << complexMatrixToString(evecInv)
	    << "Product V^{-1} * V:" << endl << tempComplexMatrixToString (evecInv_evec())
	    << "Reconstituted rate matrix:" << endl << tempComplexMatrixToString (getRateMatrix()));
}

EigenModel::~EigenModel() {
  gsl_vector_complex_free (eval);
  gsl_matrix_complex_free (evec);
  gsl_matrix_complex_free (evecInv);
}

void EigenModel::compute_exp_ev_t (double t) {
  for (AlphTok i = 0; i < model.alphabetSize(); ++i) {
    ev_t[i] = gsl_complex_mul_real (ev[i], t);
    exp_ev_t[i] = gsl_complex_exp (ev_t[i]);
  }
  LogThisAt(9,"exp(eigenvalue*" << t << "):" << complexVectorToString(exp_ev_t));
}

gsl_matrix_complex* EigenModel::getRateMatrix() const {
  gsl_matrix_complex* r = gsl_matrix_complex_alloc (model.alphabetSize(), model.alphabetSize());
  for (AlphTok i = 0; i < model.alphabetSize(); ++i)
    for (AlphTok j = 0; j < model.alphabetSize(); ++j) {
      gsl_complex rij = gsl_complex_rect (0, 0);
      for (AlphTok k = 0; k < model.alphabetSize(); ++k)
	rij = gsl_complex_add
	  (rij,
	   gsl_complex_mul (gsl_complex_mul (gsl_matrix_complex_get (evec, i, k),
					     gsl_matrix_complex_get (evecInv, k, j)),
			    ev[k]));
      gsl_matrix_complex_set (r, i, j, rij);
    }
  return r;
}

gsl_matrix_complex* EigenModel::evecInv_evec() const {
  gsl_matrix_complex* e = gsl_matrix_complex_alloc (model.alphabetSize(), model.alphabetSize());
  for (AlphTok i = 0; i < model.alphabetSize(); ++i)
    for (AlphTok j = 0; j < model.alphabetSize(); ++j) {
      gsl_complex eij = gsl_complex_rect (0, 0);
      for (AlphTok k = 0; k < model.alphabetSize(); ++k)
	eij = gsl_complex_add
	  (eij,
	   gsl_complex_mul (gsl_matrix_complex_get (evec, i, k),
			    gsl_matrix_complex_get (evecInv, k, j)));
      gsl_matrix_complex_set (e, i, j, eij);
    }
  return e;
}

double EigenModel::getSubProb (double t, AlphTok i, AlphTok j) const {
  ((EigenModel&) *this).compute_exp_ev_t (t);
  return getSubProbInner (t, i, j);
}

double EigenModel::getSubProbInner (double t, AlphTok i, AlphTok j) const {
  gsl_complex p = gsl_complex_rect (0, 0);
  for (AlphTok k = 0; k < model.alphabetSize(); ++k)
    p = gsl_complex_add
      (p,
       gsl_complex_mul (gsl_complex_mul (gsl_matrix_complex_get (evec, i, k),
					 gsl_matrix_complex_get (evecInv, k, j)),
			exp_ev_t[k]));
  Assert (SUMPROD_NEAR_REAL(p), "Probability has imaginary part: p=(%g,%g)", GSL_REAL(p), GSL_IMAG(p));
  return min (1., max (0., GSL_REAL(p)));
}

gsl_matrix* EigenModel::getSubProbMatrix (double t) const {
  gsl_matrix* sub = gsl_matrix_alloc (model.alphabetSize(), model.alphabetSize());
  ((EigenModel&) *this).compute_exp_ev_t (t);
  for (AlphTok i = 0; i < model.alphabetSize(); ++i)
    for (AlphTok j = 0; j < model.alphabetSize(); ++j)
      gsl_matrix_set (sub, i, j, getSubProbInner (t, i, j));
  return sub;
}

double EigenModel::getSubCount (AlphTok a, AlphTok b, AlphTok i, AlphTok j, const gsl_matrix* sub, const gsl_matrix_complex* eSubCount) const {
  const double p_ab = gsl_matrix_get (sub, a, b);
  const double r_ij = gsl_matrix_get (model.subRate, i, j);
  gsl_complex c_ij = gsl_complex_rect (0, 0);
  for (AlphTok k = 0; k < model.alphabetSize(); ++k) {
    gsl_complex c_ijk = gsl_complex_rect (0, 0);
    for (AlphTok l = 0; l < model.alphabetSize(); ++l)
      c_ijk = gsl_complex_add
	(c_ijk,
	 gsl_complex_mul
	 (gsl_complex_mul
	  (gsl_matrix_complex_get (evec, j, l),
	   gsl_matrix_complex_get (evecInv, l, b)),
	  gsl_matrix_complex_get (eSubCount, k, l)));
    c_ij = gsl_complex_add (c_ij,
			    gsl_complex_mul
			    (gsl_complex_mul
			     (gsl_matrix_complex_get (evec, a, k),
			      gsl_matrix_complex_get (evecInv, k, i)),
			     c_ijk));
  }
  Assert (SUMPROD_NEAR_REAL(c_ij), "Count has imaginary part: c=(%g,%g)", GSL_REAL(c_ij), GSL_IMAG(c_ij));
  return max (0., (i == j ? 1. : r_ij) * GSL_REAL(c_ij) / p_ab);
}

void EigenModel::accumSubCounts (vguard<vguard<double> >& count, AlphTok a, AlphTok b, double weight, const gsl_matrix* sub, const gsl_matrix_complex* eSubCount) const {
  for (AlphTok i = 0; i < model.alphabetSize(); ++i)
    for (AlphTok j = 0; j < model.alphabetSize(); ++j)
      count[i][j] += getSubCount (a, b, i, j, sub, eSubCount) * weight;
}

gsl_matrix_complex* EigenModel::eigenSubCount (double t) const {
  gsl_matrix_complex* esub = gsl_matrix_complex_alloc (model.alphabetSize(), model.alphabetSize());
  ((EigenModel&) *this).compute_exp_ev_t (t);
  for (AlphTok i = 0; i < model.alphabetSize(); ++i)
    for (AlphTok j = 0; j < model.alphabetSize(); ++j) {
      const bool ev_eq = i == j || SUMPROD_NEAR_EQ_COMPLEX (ev[i], ev[j]);
      gsl_matrix_complex_set
	(esub, i, j,
	 ev_eq
	 ? gsl_complex_mul_real (exp_ev_t[i], t)
	 : gsl_complex_div (gsl_complex_sub (exp_ev_t[i], exp_ev_t[j]),
			    gsl_complex_sub (ev[i], ev[j])));
    }

  LogThisAt(8,endl << "Eigensubstitution matrix at time t=" << t << ":" << endl << complexMatrixToString(esub));

  return esub;
}

vguard<vguard<double> > EigenModel::getSubCounts (const vguard<vguard<gsl_complex> >& eigenCounts) const {
  LogThisAt(8,"Eigencounts matrix:" << endl << complexMatrixToString(eigenCounts) << endl);
  vguard<vguard<double> > counts (model.alphabetSize(), vguard<double> (model.alphabetSize(), 0));
  for (AlphTok i = 0; i < model.alphabetSize(); ++i)
    for (AlphTok j = 0; j < model.alphabetSize(); ++j) {
      gsl_complex c = gsl_complex_rect (0, 0);
      for (AlphTok k = 0; k < model.alphabetSize(); ++k) {
	gsl_complex ck = gsl_complex_rect (0, 0);
	for (AlphTok l = 0; l < model.alphabetSize(); ++l)
	  ck = gsl_complex_add
	    (ck,
	     gsl_complex_mul (eigenCounts[k][l],
			      gsl_matrix_complex_get (evec, j, l)));
	c = gsl_complex_add
	  (c,
	   gsl_complex_mul (gsl_matrix_complex_get (evecInv, k, i),
			    ck));
      }
      counts[i][j] = GSL_REAL(c) * (i == j ? 1 : gsl_matrix_get (model.subRate, i, j));
    }
  return counts;
}

SumProductStorage::SumProductStorage (size_t nodes, size_t alphabetSize)
  : gappedCol (nodes),
    E (nodes, vguard<double> (alphabetSize)),
    F (nodes, vguard<double> (alphabetSize)),
    G (nodes, vguard<double> (alphabetSize)),
    logE (nodes),
    logF (nodes),
    logG (nodes)
{ }

SumProduct::SumProduct (const RateModel& model, const Tree& tree)
  : SumProductStorage (tree.nodes(), model.alphabetSize()),
    model (model),
    tree (tree),
    preorder (tree.preorderSort()),
    postorder (tree.postorderSort()),
    eigen (model),
    insProb (model.alphabetSize()),
    branchSubProb (tree.nodes(), vguard<vguard<double> > (model.alphabetSize(), vguard<double> (model.alphabetSize()))),
    branchEigenSubCount (tree.nodes())
{
  for (AlphTok i = 0; i < model.alphabetSize(); ++i)
    insProb[i] = gsl_vector_get (model.insProb, i);
  for (AlignRowIndex r = 0; r < tree.nodes() - 1; ++r) {
    ProbModel pm (model, tree.branchLength(r));
    for (AlphTok i = 0; i < model.alphabetSize(); ++i)
      for (AlphTok j = 0; j < model.alphabetSize(); ++j)
	branchSubProb[r][i][j] = gsl_matrix_get (pm.subMat, i, j);
  }

  for (AlignRowIndex r = 0; r < tree.nodes() - 1; ++r)
    branchEigenSubCount[r] = eigen.eigenSubCount (tree.branchLength(r));
}

SumProduct::~SumProduct() {
  for (auto m : branchEigenSubCount)
    gsl_matrix_complex_free (m);
}

void SumProduct::initColumn (const map<AlignRowIndex,char>& seq) {
  ungappedRows.clear();
  gappedCol = vguard<char> (tree.nodes(), Alignment::gapChar);
  vguard<int> ungappedKids (tree.nodes(), 0);
  roots.clear();
  map<size_t,SeqIdx> pos;
  for (TreeNodeIndex r = 0; r < tree.nodes(); ++r)
    if (seq.find(r) != seq.end()) {
      gappedCol[r] = seq.at(r);
      ungappedRows.push_back(r);
    }

  LogThisAt(7,"Column " << join(gappedCol,"") << " ungappedRows=(" << to_string_join(ungappedRows) << ")" << endl);
  
  for (TreeNodeIndex r = 0; r < tree.nodes(); ++r)
    if (isGap(r)) {
      fill (E[r].begin(), E[r].end(), 1);
      logE[r] = 0;
    } else {
      Require (isWild(r) || ungappedKids[r] == 0, "At node %u (%s), char %c: internal node sequences must be wildcards (%c)", r, tree.seqName(r).c_str(), seq.at(r), Alignment::wildcardChar);
      const TreeNodeIndex rp = tree.parentNode(r);
      if (rp < 0 || isGap(rp))
	roots.push_back (r);
      else
	++ungappedKids[rp];
    }
}

AlignRowIndex SumProduct::columnRoot() const {
  assertSingleRoot();
  return roots.front();
}

void SumProduct::assertSingleRoot() const {
  Require (roots.size(), "No root node in column %s tree %s", join(gappedCol,"").c_str(), tree.toString().c_str());
  Require (roots.size() == 1, "Multiple root nodes (%s) in column %s tree %s", to_string_join(roots,",").c_str(), join(gappedCol,"").c_str(), tree.toString().c_str());
}
  
void SumProduct::fillUp() {
  LogThisAt(8,"Sending tip-to-root messages, column " << join(gappedCol,"") << endl);
  colLogLike = 0;
  for (auto r : postorder) {
    logF[r] = 0;
    for (size_t nc = 0; nc < tree.nChildren(r); ++nc)
      logF[r] += logE[tree.getChild(r,nc)];
    if (!isGap(r)) {
      if (isWild(r)) {
	double Fmax = 0;
	for (AlphTok i = 0; i < model.alphabetSize(); ++i) {
	  double Fi = 1;
	  for (size_t nc = 0; nc < tree.nChildren(r); ++nc)
	    Fi *= E[tree.getChild(r,nc)][i];
	  F[r][i] = Fi;
	  if (Fi > Fmax)
	    Fmax = Fi;
	}
	if (Fmax < SUMPROD_RESCALE_THRESHOLD) {
	  for (auto& Fi: F[r])
	    Fi /= Fmax;
	  logF[r] += log (Fmax);
	}
      } else {
	for (AlphTok i = 0; i < model.alphabetSize(); ++i)
	  F[r][i] = 0;
	F[r][model.tokenize(gappedCol[r])] = 1;
      }

      const TreeNodeIndex rp = tree.parentNode(r);
      if (rp < 0 || isGap(rp))
	colLogLike += logF[r] + log (inner_product (F[r].begin(), F[r].end(), insProb.begin(), 0.));
      else {
	logE[r] = logF[r];
	for (AlphTok i = 0; i < model.alphabetSize(); ++i) {
	  double Ei = 0;
	  for (AlphTok j = 0; j < model.alphabetSize(); ++j)
	    Ei += branchSubProb[r][i][j] * F[r][j];
	  E[r][i] = Ei;
	}
      }
    }
  }
}

void SumProduct::fillDown() {
  LogThisAt(8,"Sending root-to-tip messages, column " << join(gappedCol,"") << endl);
  if (!columnEmpty()) {
    for (auto r: preorder)
      if (!isGap(r)) {
	const TreeNodeIndex rp = tree.parentNode(r);
	if (rp < 0 || isGap(rp)) {
	  G[r] = insProb;
	  logG[r] = 0;
	} else {
	  const vguard<TreeNodeIndex> rsibs = tree.getSiblings(r);
	  logG[r] = logG[rp];
	  for (auto rs: rsibs)
	    logG[r] += logE[rs];
	  for (AlphTok j = 0; j < model.alphabetSize(); ++j) {
	    double Gj = 0;
	    for (AlphTok i = 0; i < model.alphabetSize(); ++i) {
	      double p = G[rp][i] * branchSubProb[r][i][j];
	      for (auto rs: rsibs)
		if (!isGap(rs))
		  p *= E[rs][i];
	      Gj += p;
	    }
	    G[r][j] = Gj;
	  }
	}
      }
  }
}

vguard<LogProb> SumProduct::logNodePostProb (AlignRowIndex node) const {
  assertSingleRoot();
  vguard<LogProb> lpp (model.alphabetSize());
  for (AlphTok i = 0; i < model.alphabetSize(); ++i)
    lpp[i] = logF[node] + log(F[node][i]) + logG[node] + log(G[node][i]) - colLogLike;
  return lpp;
}

vguard<LogProb> SumProduct::logNodeExcludedPostProb (TreeNodeIndex node, TreeNodeIndex exclude) const {
  if (!isWild(node))
    return log_vector (F[node]);
  vguard<LogProb> lpp (model.alphabetSize(), 0);
  for (size_t nc = 0; nc < tree.nChildren(node); ++nc) {
    const TreeNodeIndex child = tree.getChild(node,nc);
    if (child != exclude)
      for (AlphTok i = 0; i < model.alphabetSize(); ++i)
	lpp[i] += log (E[child][i]);
  }
  LogProb norm = -numeric_limits<double>::infinity();
  const TreeNodeIndex parent = tree.parentNode (node);
  for (AlphTok i = 0; i < model.alphabetSize(); ++i) {
    lpp[i] += log (parent == exclude ? insProb[i] : G[node][i]);
    log_accum_exp (norm, lpp[i]);
  }
  for (auto& lp: lpp)
    lp -= norm;
  return lpp;
}

LogProb SumProduct::logBranchPostProb (AlignRowIndex node, AlphTok parentState, AlphTok nodeState) const {
  assertSingleRoot();
  const TreeNodeIndex parent = tree.parentNode(node);
  const TreeNodeIndex sibling = tree.getSibling(node);
  return logG[parent] + log(G[parent][parentState]) + log(branchSubProb[node][parentState][nodeState]) + logF[node] + log(F[node][nodeState]) + logE[sibling] + log(E[sibling][parentState]) - colLogLike;
}

AlphTok SumProduct::maxPostState (AlignRowIndex node) const {
  const auto lpp = logNodePostProb (node);
  return (AlphTok) (max_element(lpp.begin(),lpp.end()) - lpp.begin());
}

void SumProduct::accumulateRootCounts (vguard<double>& rootCounts, double weight) const {
  const auto rootNode = columnRoot();
  const double norm = exp (logF[rootNode] - colLogLike);
  for (AlphTok i = 0; i < model.alphabetSize(); ++i)
    rootCounts[i] += weight * insProb[i] * F[rootNode][i] * norm;
}

void SumProduct::accumulateSubCounts (vguard<double>& rootCounts, vguard<vguard<double> >& subCounts, double weight) const {
  LogThisAt(8,"Accumulating substitution counts, column " << join(gappedCol,"") << endl);
  accumulateRootCounts (rootCounts, weight);

  const auto rootNode = columnRoot();
  for (auto node : ungappedRows)
    if (node != rootNode) {
      LogThisAt(9,"Accumulating substitution counts, column " << join(gappedCol,"") << " node " << tree.seqName(node) << endl);
      const TreeNodeIndex parent = tree.parentNode(node);
      gsl_matrix* submat = eigen.getSubProbMatrix (tree.branchLength(node));
      for (AlphTok a = 0; a < model.alphabetSize(); ++a)
	for (AlphTok b = 0; b < model.alphabetSize(); ++b)
	  eigen.accumSubCounts (subCounts, a, b, weight * exp (logBranchPostProb (node, a, b)), submat, branchEigenSubCount[node]);
      gsl_matrix_free (submat);
    }
}

void SumProduct::accumulateEigenCounts (vguard<double>& rootCounts, vguard<vguard<gsl_complex> >& eigenCounts, double weight) const {
  LogThisAt(8,"Accumulating eigencounts, column " << join(gappedCol,"") << endl);
  accumulateRootCounts (rootCounts, weight);

  const auto rootNode = columnRoot();
  const int A = model.alphabetSize();
  vguard<double> U (A), D (A);
  vguard<gsl_complex> Ubasis (A), Dbasis (A);
  vguard<double> U0 (A), D0 (A);
  for (auto node : ungappedRows)
    if (node != rootNode) {
      LogThisAt(9,"Accumulating eigencounts, column " << join(gappedCol,"") << " node " << tree.seqName(node) << endl);
      const TreeNodeIndex parent = tree.parentNode(node);
      const TreeNodeIndex sibling = tree.getSibling(node);
      const vguard<double>& U0 = F[node];
      for (AlphTok i = 0; i < A; ++i)
	D0[i] = G[parent][i] * E[sibling][i];
      const double maxU0 = *max_element (U0.begin(), U0.end());
      const double maxD0 = *max_element (D0.begin(), D0.end());
      const double norm = exp (colLogLike - logF[node] - logG[parent] - logE[sibling]) / (maxU0 * maxD0);

      // U[b] = U0[b] / maxU0; Ubasis[l] = sum_b U[b] * evecInv[l][b]
      for (AlphTok b = 0; b < A; ++b)
	U[b] = U0[b] / maxU0;

      for (AlphTok l = 0; l < A; ++l) {
	Ubasis[l] = gsl_complex_rect (0, 0);
	for (AlphTok b = 0; b < A; ++b)
	  Ubasis[l] = gsl_complex_add
	    (Ubasis[l],
	     gsl_complex_mul_real (gsl_matrix_complex_get (eigen.evecInv, l, b),
				   U[b]));
      }

      // D[a] = D0[a] / maxD0; Dbasis[k] = sum_a D[a] * evec[a][k]
      for (AlphTok a = 0; a < A; ++a)
	D[a] = D0[a] / maxD0;

      for (AlphTok k = 0; k < A; ++k) {
	Dbasis[k] = gsl_complex_rect (0, 0);
	for (AlphTok a = 0; a < A; ++a)
	  Dbasis[k] = gsl_complex_add
	    (Dbasis[k],
	     gsl_complex_mul_real (gsl_matrix_complex_get (eigen.evec, a, k),
				   D[a]));
      }

      // R = evec * evals * evecInv
      // exp(RT) = evec * exp(evals T) * evecInv
      // count(i,j|a,b,T) = Q / exp(RT)_ab
      // where...
      // Q = \sum_a \sum_b \int_{t=0}^T D_a exp(Rt)_ai R_ij exp(R(T-t))_jb U_b dt
      //   = \sum_a \sum_b \int_{t=0}^T D_a (\sum_k evec_ak exp(eval_k t) evecInv_ki) R_ij (\sum_l evec_jl exp(eval_l (T-t)) evecInv_lb) U_b dt
      //   = \sum_a \sum_b D_a \sum_k evec_ak evecInv_ki R_ij \sum_l evec_jl evecInv_lb U_b \int_{t=0}^T exp(eval_k t) exp(eval_l (T-t)) dt
      //   = \sum_a \sum_b D_a \sum_k evec_ak evecInv_ki R_ij \sum_l evec_jl evecInv_lb U_b eigenSubCount(k,l,T)
      //   = R_ij \sum_k evecInv_ki \sum_l evec_jl (\sum_a D_a evec_ak) (\sum_b U_b evecInv_lb) eigenSubCount(k,l,T)
      //   = R_ij \sum_k evecInv_ki \sum_l evec_jl Dbasis_k Ubasis_l eigenSubCount(k,l,T)

      // eigenCounts[k][l] += Dbasis[k] * eigenSub[k][l] * Ubasis[l] / norm
      for (AlphTok k = 0; k < A; ++k)
	for (AlphTok l = 0; l < A; ++l)
	  eigenCounts[k][l] =
	    gsl_complex_add
	    (eigenCounts[k][l],
	     gsl_complex_mul_real
	     (gsl_complex_mul
	      (Dbasis[k],
	       gsl_complex_mul
	       (gsl_matrix_complex_get (branchEigenSubCount[node], k, l),
		Ubasis[l])),
	      weight / norm));
    }
}

string tempComplexMatrixToString (gsl_matrix_complex* mx) {
  const string s = complexMatrixToString (mx);
  gsl_matrix_complex_free (mx);
  return s;
}

string complexMatrixToString (const gsl_matrix_complex* mx) {
  ostringstream s;
  for (size_t i = 0; i < mx->size1; ++i) {
    for (size_t j = 0; j < mx->size2; ++j) {
      const gsl_complex c = gsl_matrix_complex_get (mx, i, j);
      s << " (" << GSL_REAL(c) << "," << GSL_IMAG(c) << ")";
    }
    s << endl;
  }
  return s.str();
}

string complexMatrixToString (const vguard<vguard<gsl_complex> >& mx) {
  ostringstream s;
  for (size_t i = 0; i < mx.size(); ++i) {
    for (size_t j = 0; j < mx[i].size(); ++j) {
      const gsl_complex& c = mx[i][j];
      s << " (" << GSL_REAL(c) << "," << GSL_IMAG(c) << ")";
    }
    s << endl;
  }
  return s.str();
}

string complexVectorToString (const gsl_vector_complex* v) {
  ostringstream s;
  for (size_t i = 0; i < v->size; ++i) {
    const gsl_complex c = gsl_vector_complex_get (v, i);
      s << " (" << GSL_REAL(c) << "," << GSL_IMAG(c) << ")";
  }
  s << endl;
  return s.str();
}

string complexVectorToString (const vector<gsl_complex>& v) {
  ostringstream s;
  for (size_t i = 0; i < v.size(); ++i) {
    const gsl_complex c = v[i];
      s << " (" << GSL_REAL(c) << "," << GSL_IMAG(c) << ")";
  }
  s << endl;
  return s.str();
}

AlignColSumProduct::AlignColSumProduct (const RateModel& model, const Tree& tree, const vguard<FastSeq>& gapped)
  : SumProduct (model, tree),
    gapped (gapped),
    col (0)
{
  Require (tree.nodes() == gapped.size(), "Every tree node must have an alignment row");
  initAlignColumn();
}

void AlignColSumProduct::initAlignColumn() {
  map<AlignRowIndex,char> seq;
  for (TreeNodeIndex r = 0; r < tree.nodes(); ++r)
    if (!Alignment::isGap (gapped[r].seq[col]))
      seq[r] = gapped[r].seq[col];
  initColumn (seq);
}

bool AlignColSumProduct::alignmentDone() const {
  return col >= gapped.front().length();
}

void AlignColSumProduct::nextColumn() {
  ++col;
  if (!alignmentDone())
    initAlignColumn();
}

void AlignColSumProduct::appendAncestralReconstructedColumn (vguard<FastSeq>& out) const {
  if (col == 0) {
    out = gapped;
    for (auto& fs : out) {
      fs.seq.clear();
      fs.qual.clear();
    }
  }
  for (AlignRowIndex row = 0; row < gapped.size(); ++row) {
    const char g = gapped[row].seq[col];
    out[row].seq.push_back (Alignment::isWildcard(g) ? model.alphabet[maxPostState(row)] : g);
  }
}

void AlignColSumProduct::appendAncestralPostProbColumn (ReconPostProbMap& rpp, double minProb, double maxProb) const {
  const LogProb lpMin = log(minProb), lpMax = log(maxProb);
  for (AlignRowIndex row = 0; row < gapped.size(); ++row) {
    const char g = gapped[row].seq[col];
    if (Alignment::isWildcard(g)) {
      auto lp = logNodePostProb (row);
      for (AlphTok tok = 0; tok < model.alphabet.size(); ++tok)
	if (lp[tok] >= lpMin && lp[tok] <= lpMax)
	  rpp[row][col][model.alphabet[tok]] = exp(lp[tok]);
    }
  }
}
