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
  transform (alphabet.begin(), alphabet.end(), alphabet.begin(), ::toupper);
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
  gsl_matrix* QR = newAlphabetMatrix();
  gsl_vector* tau = newAlphabetVector();
  gsl_vector* zero = newAlphabetVector();
  gsl_vector* eqm = newAlphabetVector();

  gsl_vector_set_zero (zero);
  gsl_vector_set_zero (eqm);

  CheckGsl (gsl_matrix_memcpy (QR, subRate));
  CheckGsl (gsl_linalg_QR_decomp (QR, tau));

  fprintf (stderr, "QR:\n");
  gsl_matrix_fprintf (stderr, QR, "%g");
  fprintf (stderr, "tau:\n");
  gsl_vector_fprintf (stderr, tau, "%g");
  
  CheckGsl (gsl_linalg_QR_solve (QR, tau, zero, eqm));

  fprintf (stderr, "eqm:\n");
  gsl_vector_fprintf (stderr, eqm, "%g");

  double eqmNorm = 0;
  for (AlphTok i = 0; i < alphabetSize(); ++i)
    eqmNorm += gsl_vector_get (eqm, i);
  Assert (eqmNorm != 0, "Equilibrium vector is zero");

  CheckGsl (gsl_vector_scale (eqm, 1 / eqmNorm));

  gsl_vector_free (tau);
  gsl_vector_free (zero);
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

