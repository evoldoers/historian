#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_math.h>

#include "sumprod.h"
#include "util.h"
#include "logger.h"

#define SUMPROD_RESCALE_THRESHOLD 1e-30

SumProductStorage::SumProductStorage (size_t components, size_t nodes, size_t alphabetSize)
  : gappedCol (nodes),
    E (components, vguard<vguard<double> > (nodes, vguard<double> (alphabetSize))),
    F (components, vguard<vguard<double> > (nodes, vguard<double> (alphabetSize))),
    G (components, vguard<vguard<double> > (nodes, vguard<double> (alphabetSize))),
    logE (components, vguard<double> (nodes)),
    logF (components, vguard<double> (nodes)),
    logG (components, vguard<double> (nodes)),
    cptLogLike (components)
{ }

SumProduct::SumProduct (const RateModel& model, const Tree& tree)
  : SumProductStorage (model.components(), tree.nodes(), model.alphabetSize()),
    model (model),
    tree (tree),
    preorder (tree.preorderSort()),
    postorder (tree.postorderSort()),
    eigen (model),
    insProb (model.components(), vguard<double> (model.alphabetSize())),
    branchSubProb (model.components(), vguard<vguard<vguard<double> > > (tree.nodes(), vguard<vguard<double> > (model.alphabetSize(), vguard<double> (model.alphabetSize())))),
    branchEigenSubCount (model.components(), vguard<gsl_matrix_complex*> (tree.nodes())),
    logCptWeight (log_vector (model.cptWeight))
{
  for (int cpt = 0; cpt < components(); ++cpt)
    for (AlphTok i = 0; i < model.alphabetSize(); ++i)
      insProb[cpt][i] = gsl_vector_get (model.insProb[cpt], i);

  for (AlignRowIndex r = 0; r < tree.nodes() - 1; ++r) {
      ProbModel pm (model, tree.branchLength(r));
      for (int cpt = 0; cpt < components(); ++cpt)
	for (AlphTok i = 0; i < model.alphabetSize(); ++i)
	  for (AlphTok j = 0; j < model.alphabetSize(); ++j)
	    branchSubProb[cpt][r][i][j] = gsl_matrix_get (pm.subMat[cpt], i, j);
  }

  for (AlignRowIndex r = 0; r < tree.nodes() - 1; ++r) {
    vguard<gsl_matrix_complex*> esc = eigen.eigenSubCount (tree.branchLength(r));
    for (int cpt = 0; cpt < components(); ++cpt)
      branchEigenSubCount[cpt][r] = esc[cpt];
  }
}

SumProduct::~SumProduct() {
  for (int cpt = 0; cpt < components(); ++cpt)
    for (auto m : branchEigenSubCount[cpt])
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
      const char c = seq.at(r);
      gappedCol[r] = model.isValidSymbol(c) ? c : Alignment::wildcardChar;
      ungappedRows.push_back(r);
    }

  LogThisAt(7,"Column " << join(gappedCol,"") << " ungappedRows=(" << to_string_join(ungappedRows) << ")" << endl);
  
  for (TreeNodeIndex r = 0; r < tree.nodes(); ++r)
    if (isGap(r)) {
      for (int cpt = 0; cpt < components(); ++cpt) {
	fill (E[cpt][r].begin(), E[cpt][r].end(), 1);
	logE[cpt][r] = 0;
      }
    } else {
      //      Require (isWild(r) || ungappedKids[r] == 0, "At node %u (%s), char %c: internal node sequences must be wildcards (%c)", r, tree.seqName(r).c_str(), seq.at(r), Alignment::wildcardChar);
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
  colLogLike = -numeric_limits<double>::infinity();
  for (int cpt = 0; cpt < components(); ++cpt) {
    LogThisAt(8,"Sending tip-to-root messages, component #" << cpt << " column " << join(gappedCol,"") << endl);
    cptLogLike[cpt] = 0;
    for (auto r : postorder) {
      logF[cpt][r] = 0;
      for (size_t nc = 0; nc < tree.nChildren(r); ++nc)
	logF[cpt][r] += logE[cpt][tree.getChild(r,nc)];
      if (!isGap(r)) {
	const char c = gappedCol[r];
	if (Alignment::isWildcard(c)) {
	  double Fmax = 0;
	  for (AlphTok i = 0; i < model.alphabetSize(); ++i) {
	    double Fi = 1;
	    for (size_t nc = 0; nc < tree.nChildren(r); ++nc)
	      Fi *= E[cpt][tree.getChild(r,nc)][i];
	    F[cpt][r][i] = Fi;
	    if (Fi > Fmax)
	      Fmax = Fi;
	  }
	  if (Fmax < SUMPROD_RESCALE_THRESHOLD) {
	    for (auto& Fi: F[cpt][r])
	      Fi /= Fmax;
	    logF[cpt][r] += log (Fmax);
	  }
	} else {  // !isWild(r)
	  const AlphTok tok = model.tokenize(c);
	  double Ftok = 1;
	  for (size_t nc = 0; nc < tree.nChildren(r); ++nc)
	    Ftok *= E[cpt][tree.getChild(r,nc)][tok];

	  if (Ftok < SUMPROD_RESCALE_THRESHOLD) {
	    logF[cpt][r] += log (Ftok);
	    Ftok = 1;
	  }

	  for (AlphTok i = 0; i < model.alphabetSize(); ++i)
	    F[cpt][r][i] = 0;
	  F[cpt][r][tok] = Ftok;
	}

	LogThisAt(10,"Row " << setw(3) << r << " " << c << " " << cpt
		  << " logF=" << setw(9) << setprecision(3) << logF[cpt][r]
		  << " F=(" << to_string_join(F[cpt][r]," ",9,3) << ")" << endl);

	const TreeNodeIndex rp = tree.parentNode(r);
	if (rp < 0 || isGap(rp))
	  cptLogLike[cpt] += logF[cpt][r] + log (inner_product (F[cpt][r].begin(), F[cpt][r].end(), insProb[cpt].begin(), 0.));
	else {
	  logE[cpt][r] = logF[cpt][r];
	  for (AlphTok i = 0; i < model.alphabetSize(); ++i) {
	    double Ei = 0;
	    for (AlphTok j = 0; j < model.alphabetSize(); ++j)
	      Ei += branchSubProb[cpt][r][i][j] * F[cpt][r][j];
	    E[cpt][r][i] = Ei;
	  }
	}
      }
    }
    log_accum_exp (colLogLike, logCptWeight[cpt] + cptLogLike[cpt]);
  }
}

void SumProduct::fillDown() {
  for (int cpt = 0; cpt < components(); ++cpt) {
    LogThisAt(8,"Sending root-to-tip messages, component #" << cpt << " column " << join(gappedCol,"") << endl);
    if (!columnEmpty()) {
      for (auto r: preorder) {
	if (!isGap(r)) {
	  const TreeNodeIndex rp = tree.parentNode(r);
	  if (rp < 0 || isGap(rp)) {
	    G[cpt][r] = insProb[cpt];
	    logG[cpt][r] = 0;
	  } else {
	    const vguard<TreeNodeIndex> rsibs = tree.getSiblings(r);
	    logG[cpt][r] = logG[cpt][rp];
	    for (auto rs: rsibs)
	      logG[cpt][r] += logE[cpt][rs];
	    for (AlphTok j = 0; j < model.alphabetSize(); ++j) {
	      double Gj = 0;
	      for (AlphTok i = 0; i < model.alphabetSize(); ++i) {
		double p = G[cpt][rp][i] * branchSubProb[cpt][r][i][j];
		for (auto rs: rsibs)
		  if (!isGap(rs))
		    p *= E[cpt][rs][i];
		Gj += p;
	      }
	      G[cpt][r][j] = Gj;
	    }
	  }
	}

	LogThisAt(10,"Row " << setw(3) << r << " " << cpt << " " << gappedCol[r]
		  << " logG=" << setw(9) << setprecision(3) << logG[cpt][r]
		  << " G=(" << to_string_join(G[cpt][r]," ",9,3) << ")" << endl);
      }
    }
  }
}

LogProb SumProduct::computeColumnLogLikelihoodAt (AlignRowIndex node) const {
  LogProb lp = -numeric_limits<double>::infinity();
  for (int cpt = 0; cpt < components(); ++cpt)
    for (AlphTok i = 0; i < model.alphabetSize(); ++i)
      log_accum_exp (lp, logCptWeight[cpt] + logF[cpt][node] + log(F[cpt][node][i]) + logG[cpt][node] + log(G[cpt][node][i]));
  return lp;
}

vguard<LogProb> SumProduct::logNodePostProb (AlignRowIndex node) const {
  assertSingleRoot();
  vguard<LogProb> lpp (model.alphabetSize(), -numeric_limits<double>::infinity());
  for (AlphTok i = 0; i < model.alphabetSize(); ++i) {
    for (int cpt = 0; cpt < components(); ++cpt)
      log_accum_exp (lpp[i], logCptWeight[cpt] + logF[cpt][node] + log(F[cpt][node][i]) + logG[cpt][node] + log(G[cpt][node][i]) - colLogLike);
    lpp[i] = min (lpp[i], 0.);  // guard against overflow probability > 1
  }
  return lpp;
}

vguard<vguard<LogProb> > SumProduct::logNodeExcludedPostProb (TreeNodeIndex node, TreeNodeIndex exclude, bool normalize) const {
  Require (!isGap(node), "Attempt to find posterior probability of sequence at gapped position");
  const UnvalidatedAlphTok tok = isWild(node) ? -1 : model.tokenize(gappedCol[node]);
  vguard<LogProb> lppInit (model.alphabetSize(), isWild(node) ? 0 : -numeric_limits<double>::infinity());
  if (!isWild(node))
    lppInit[tok] = 0;
  vguard<vguard<LogProb> > v (model.components(), lppInit);
  LogProb norm = -numeric_limits<double>::infinity();
  for (int cpt = 0; cpt < components(); ++cpt) {
    vguard<LogProb>& lpp = v[cpt];
    for (auto& lp: lpp)
      lp += logCptWeight[cpt];
    for (size_t nc = 0; nc < tree.nChildren(node); ++nc) {
      const TreeNodeIndex child = tree.getChild(node,nc);
      if (child != exclude)
	for (AlphTok i = 0; i < model.alphabetSize(); ++i)
	  lpp[i] += log (E[cpt][child][i]) + logE[cpt][child];
    }
    const TreeNodeIndex parent = tree.parentNode (node);
    for (AlphTok i = 0; i < model.alphabetSize(); ++i) {
      lpp[i] += parent == exclude
	? 0 // to add a prior for orphaned nodes, this should be log(insProb[i]), but that complicates MCMC etc
	: (log(G[cpt][node][i]) + logG[cpt][node]);
      log_accum_exp (norm, lpp[i]);
    }
  }
  if (normalize)
    for (auto& lpp: v)
      for (auto& lp: lpp)
	lp -= norm;
  return v;
}

LogProb SumProduct::logBranchPostProb (int cpt, AlignRowIndex node, AlphTok parentState, AlphTok nodeState) const {
  assertSingleRoot();
  const TreeNodeIndex parent = tree.parentNode(node);
  const TreeNodeIndex sibling = tree.getSibling(node);
  return logCptWeight[cpt] + logG[cpt][parent] + log(G[cpt][parent][parentState]) + log(branchSubProb[cpt][node][parentState][nodeState]) + logF[cpt][node] + log(F[cpt][node][nodeState]) + logE[cpt][sibling] + log(E[cpt][sibling][parentState]) - colLogLike;
}

AlphTok SumProduct::maxPostState (AlignRowIndex node) const {
  const auto lpp = logNodePostProb (node);
  return (AlphTok) (max_element(lpp.begin(),lpp.end()) - lpp.begin());
}

void SumProduct::accumulateRootCounts (vguard<vguard<double> >& rootCounts, double weight) const {
  const auto rootNode = columnRoot();
  for (int cpt = 0; cpt < components(); ++cpt) {
    const double norm = exp (logCptWeight[cpt] + logF[cpt][rootNode] - colLogLike);
    for (AlphTok i = 0; i < model.alphabetSize(); ++i)
      rootCounts[cpt][i] += weight * insProb[cpt][i] * F[cpt][rootNode][i] * norm;
  }
}

void SumProduct::accumulateSubCounts (vguard<vguard<double> >& rootCounts, vguard<vguard<vguard<double> > >& subCounts, double weight) const {
  LogThisAt(8,"Accumulating substitution counts, column " << join(gappedCol,"") << ", weight " << weight << endl);
  accumulateRootCounts (rootCounts, weight);

  const auto rootNode = columnRoot();
  for (auto node : ungappedRows)
    if (node != rootNode) {
      LogThisAt(9,"Accumulating substitution counts, column " << join(gappedCol,"") << " node " << tree.seqName(node) << endl);
      const TreeNodeIndex parent = tree.parentNode(node);
      vguard<gsl_matrix*> submat = model.getSubProbMatrix (tree.branchLength(node));
      for (int cpt = 0; cpt < components(); ++cpt) {
	LogThisAt(9,"Accumulating substitution counts, column " << join(gappedCol,"") << " node " << tree.seqName(node) << " component #" << cpt << endl);
	for (AlphTok a = 0; a < model.alphabetSize(); ++a)
	  for (AlphTok b = 0; b < model.alphabetSize(); ++b)
	    eigen.accumSubCounts (cpt, subCounts[cpt], a, b, weight * exp (logBranchPostProb (cpt, node, a, b)), submat[cpt], branchEigenSubCount[cpt][node]);
      }
      for (auto& sm: submat)
	gsl_matrix_free (sm);
    }
}

void SumProduct::accumulateEigenCounts (vguard<vguard<double> >& rootCounts, vguard<vguard<vguard<gsl_complex> > >& eigenCounts, double weight) const {
  LogThisAt(8,"Accumulating eigencounts, column " << join(gappedCol,"") << ", weight " << weight << endl);
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
      for (int cpt = 0; cpt < components(); ++cpt) {
	LogThisAt(9,"Accumulating eigencounts, column " << join(gappedCol,"") << " node " << tree.seqName(node) << " component #" << cpt << endl);
	const vguard<double>& U0 = F[cpt][node];
	for (AlphTok i = 0; i < A; ++i)
	  D0[i] = G[cpt][parent][i] * E[cpt][sibling][i];
	const double maxU0 = *max_element (U0.begin(), U0.end());
	const double maxD0 = *max_element (D0.begin(), D0.end());
	const double norm = exp (colLogLike - logCptWeight[cpt] - logF[cpt][node] - logG[cpt][parent] - logE[cpt][sibling]) / (maxU0 * maxD0);

	// U[b] = U0[b] / maxU0; Ubasis[l] = sum_b U[b] * evecInv[l][b]
	for (AlphTok b = 0; b < A; ++b)
	  U[b] = U0[b] / maxU0;

	for (AlphTok l = 0; l < A; ++l) {
	  Ubasis[l] = gsl_complex_rect (0, 0);
	  for (AlphTok b = 0; b < A; ++b)
	    Ubasis[l] = gsl_complex_add
	      (Ubasis[l],
	       gsl_complex_mul_real (gsl_matrix_complex_get (eigen.evecInv[cpt], l, b),
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
	       gsl_complex_mul_real (gsl_matrix_complex_get (eigen.evec[cpt], a, k),
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
	    eigenCounts[cpt][k][l] =
	      gsl_complex_add
	      (eigenCounts[cpt][k][l],
	       gsl_complex_mul_real
	       (gsl_complex_mul
		(Dbasis[k],
		 gsl_complex_mul
		 (gsl_matrix_complex_get (branchEigenSubCount[cpt][node], k, l),
		  Ubasis[l])),
		weight / norm));

	//	LogThisAt(9,"colLogLike=" << colLogLike << " logF[cpt][node]=" << logF[cpt][node] << " logG[cpt][parent]=" << logG[cpt][parent] << " logE[cpt][sibling]=" << logE[cpt][sibling] << " maxU0=" << maxU0 << " maxD0=" << maxD0 << endl);
	//	LogThisAt(8,"Component #" << cpt << " eigencounts matrix (norm=" << norm << "):" << endl << complexMatrixToString(eigenCounts[cpt]) << endl);
      }
    }
}

AlignColSumProduct::AlignColSumProduct (const RateModel& model, const Tree& tree, const vguard<FastSeq>& gapped)
  : SumProduct (model, tree),
    gapped (gapped),
    col (0)
{
  Assert (tree.nodes() == gapped.size(), "Number of nodes in tree (%d) does not match number of sequences (%d)", tree.nodes(), gapped.size());
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
