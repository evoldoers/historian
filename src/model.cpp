#include <math.h>
#include <gsl/gsl_linalg.h>
#include <algorithm>
#include <set>

#include "model.h"
#include "jsonutil.h"
#include "util.h"

void AlphabetOwner::readAlphabet (const JsonValue& json) {
  const JsonMap jm (json);
  Assert (jm.containsType("alphabet",JSON_STRING), "No alphabet");
  alphabet = jm["alphabet"].toString();
  set<char> s;
  for (auto c : alphabet) {
    Assert (s.find(c) == s.end(), "Duplicate character %c in alphabet %s", c, alphabet.c_str());
    s.insert (c);
  }
}

UnvalidatedAlphTok AlphabetOwner::tokenize (char c) const {
  return ::tokenize (c, alphabet);
}

gsl_matrix* AlphabetOwner::newAlphabetMatrix() const {
  return gsl_matrix_alloc (alphabetSize(), alphabetSize());
}

gsl_vector* AlphabetOwner::newAlphabetVector() const {
  return gsl_vector_alloc (alphabetSize());
}

RateModel::RateModel()
  : subRate(NULL)
{ }

void RateModel::read (const JsonValue& json) {
  Assert (subRate == NULL, "RateModel already initialized");

  readAlphabet (json);

  subRate = newAlphabetMatrix();
  gsl_matrix_set_zero (subRate);

  const JsonMap jm (json);
  const JsonMap rateMatrix = jm.getObject ("subrate");
  for (auto ci : alphabet) {
    AlphTok i = (AlphTok) tokenize (ci);
    const string si (1, ci);
    if (rateMatrix.contains (si)) {
      const JsonMap rateMatrixRow = rateMatrix.getObject (si);
      for (auto cj : alphabet)
	if (cj != ci) {
	  AlphTok j = (AlphTok) tokenize (cj);
	  const string sj (1, cj);
	  if (rateMatrixRow.contains (sj)) {
	    const double rate = rateMatrixRow.getNumber (string (1, cj));
	    *(gsl_matrix_ptr (subRate, i, j)) += rate;
	    *(gsl_matrix_ptr (subRate, i, i)) -= rate;
	  }
	}
    }
  }

  insRate = jm.getNumber("insrate");
  insExtProb = jm.getNumber("insextprob");
  delRate = jm.getNumber("delrate");
  delExtProb = jm.getNumber("delextprob");
}

void RateModel::write (ostream& out) const {
  out << "{" << endl;
  out << " \"alphabet\": \"" << alphabet << "\"," << endl;
  out << " \"insrate\": " << insRate << "," << endl;
  out << " \"insextprob\": " << insExtProb << "," << endl;
  out << " \"delrate\": " << delRate << "," << endl;
  out << " \"delextprob\": " << delExtProb << "," << endl;
  out << " \"subrate\": {" << endl;
  for (AlphTok i = 0; i < alphabetSize(); ++i) {
    out << "  \"" << alphabet[i] << "\": {" << endl;
    for (AlphTok j = 0; j < alphabetSize(); ++j)
      if (i != j)
	out << "   \"" << alphabet[j] << "\": " << gsl_matrix_get(subRate,i,j) << (j < alphabetSize() - (i == alphabetSize() - 1 ? 2 : 1) ? ",\n" : "");
    out << endl << "  }" << (i < alphabetSize() - 1 ? "," : "") << endl;
  }
  out << " }" << endl;
  out << "}" << endl;
}

gsl_vector* RateModel::getEqmProb() const {
  // find eqm via QR decomposition
  gsl_matrix* QR = gsl_matrix_alloc (alphabetSize() + 1, alphabetSize());
  for (AlphTok j = 0; j < alphabetSize(); ++j) {
    for (AlphTok i = 0; i < alphabetSize(); ++i)
      gsl_matrix_set (QR, i, j, gsl_matrix_get (subRate, j, i));
    gsl_matrix_set (QR, alphabetSize(), j, 1);
  }

  gsl_vector* tau = newAlphabetVector();

  CheckGsl (gsl_linalg_QR_decomp (QR, tau));

  gsl_vector* b = gsl_vector_alloc (alphabetSize() + 1);
  gsl_vector_set_zero (b);
  gsl_vector_set (b, alphabetSize(), 1);

  gsl_vector* residual = gsl_vector_alloc (alphabetSize() + 1);
  gsl_vector* eqm = newAlphabetVector();
  
  CheckGsl (gsl_linalg_QR_lssolve (QR, tau, b, eqm, residual));

  // make sure no "probabilities" have become negative & (re-)normalize
  double eqmNorm = 0;
  for (AlphTok i = 0; i < alphabetSize(); ++i) {
    const double eqm_i = max (0., gsl_vector_get (eqm, i));
    gsl_vector_set (eqm, i, eqm_i);
    eqmNorm += eqm_i;
  }
  CheckGsl (gsl_vector_scale (eqm, 1 / eqmNorm));

  gsl_vector_free (tau);
  gsl_vector_free (b);
  gsl_vector_free (residual);
  gsl_matrix_free (QR);
  
  return eqm;
}

gsl_matrix* RateModel::getSubProb (double t) const {
  gsl_matrix* m = newAlphabetMatrix();
  gsl_matrix* rt = newAlphabetMatrix();
  CheckGsl (gsl_matrix_memcpy (rt, subRate));
  CheckGsl (gsl_matrix_scale (rt, t));
  CheckGsl (gsl_linalg_exponential_ss (rt, m, GSL_PREC_DOUBLE));
  gsl_matrix_free (rt);
  return m;
}

ProbModel::ProbModel (const RateModel& model, double t)
  : AlphabetOwner (model),
    t (t),
    ins (1 - exp (-model.insRate * t)),
    del (1 - exp (-model.delRate * t)),
    insExt (model.insExtProb),
    delExt (model.delExtProb),
    insVec (model.getEqmProb()),
    subMat (model.getSubProb (t))
{ }

ProbModel::~ProbModel() {
  gsl_matrix_free (subMat);
  gsl_vector_free (insVec);
}

void ProbModel::write (ostream& out) const {
  out << "{" << endl;
  out << " \"alphabet\": \"" << alphabet << "\"," << endl;
  out << " \"insBegin\": " << ins << "," << endl;
  out << " \"insExtend\": " << insExt << "," << endl;
  out << " \"delBegin\": " << del << "," << endl;
  out << " \"delExtend\": " << delExt << "," << endl;
  out << " \"insVec\": {" << endl;
  for (AlphTok i = 0; i < alphabetSize(); ++i)
    out << "  \"" << alphabet[i] << "\": " << gsl_vector_get(insVec,i) << (i < alphabetSize() - 1 ? ",\n" : "");
  out << endl << " }," << endl;
  out << " \"subMat\": {" << endl;
  for (AlphTok i = 0; i < alphabetSize(); ++i) {
    out << "  \"" << alphabet[i] << "\": {" << endl;
    for (AlphTok j = 0; j < alphabetSize(); ++j)
      out << "   \"" << alphabet[j] << "\": " << gsl_matrix_get(subMat,i,j) << (j < alphabetSize() - 1 ? ",\n" : "");
    out << endl << "  }" << (i < alphabetSize() - 1 ? "," : "") << endl;
  }
  out << " }" << endl;
  out << "}" << endl;
}

LogProbModel::LogProbModel (const ProbModel& pm)
  : logInsProb (pm.alphabetSize()),
    logSubProb (pm.alphabetSize(), vguard<LogProb> (pm.alphabetSize()))
{
  for (AlphTok i = 0; i < pm.alphabetSize(); ++i) {
    logInsProb[i] = log (gsl_vector_get (pm.insVec, i));
    for (AlphTok j = 0; j < pm.alphabetSize(); ++j)
      logSubProb[i][j] = log (gsl_matrix_get (pm.subMat, i, j));
  }
}
