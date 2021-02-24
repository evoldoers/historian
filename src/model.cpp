#include <math.h>

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_complex_math.h>

#include <iomanip>
#include <algorithm>
#include <set>

#include "model.h"
#include "jsonutil.h"
#include "util.h"
#include "alignpath.h"
#include "logger.h"
#include "sumprod.h"

#define EIGENMODEL_EPSILON 1e-6
#define EIGENMODEL_NEAR_EQ(X,Y) (gsl_fcmp (X, Y, EIGENMODEL_EPSILON) == 0)
#define EIGENMODEL_NEAR_EQ_COMPLEX(X,Y) (EIGENMODEL_NEAR_EQ(GSL_REAL(X),GSL_REAL(Y)) && EIGENMODEL_NEAR_EQ(GSL_IMAG(X),GSL_IMAG(Y)))
#define EIGENMODEL_NEAR_REAL(X) (abs(GSL_IMAG(X)) < EIGENMODEL_EPSILON)

struct DistanceMatrixParams {
  const map<pair<AlphTok,AlphTok>,int>& pairCount;
  const RateModel& model;
  DistanceMatrixParams (const map<pair<AlphTok,AlphTok>,int>& pairCount, const RateModel& model)
    : pairCount(pairCount),
      model(model)
  { }
  double tJC() const;  // Jukes-Cantor estimate
  double tML (int maxIterations) const;
  double negLogLike (double t) const;
};

void AlphabetOwner::readAlphabet (const JsonValue& json) {
  const JsonMap jm (json);
  Assert (jm.containsType("alphabet",JSON_STRING), "No alphabet");
  initAlphabet (jm["alphabet"].toString(),
		jm.containsType("wildcard",JSON_STRING) ? jm["wildcard"].toString()[0] : Alignment::wildcardChar);
}

void AlphabetOwner::initAlphabet (const string& a, char wild) {
  alphabet = a;
  wildcard = wild;
  set<char> s;
  for (auto c : alphabet) {
    Assert (s.find(c) == s.end(), "Duplicate character %c in alphabet %s", c, alphabet.c_str());
    Assert (c != Alignment::wildcardChar, "Character %c is reserved for internal use as a wildcard, and cannot be used as a character in alphabet %s", c, alphabet.c_str());
    Assert (c != Alignment::gapChar, "Character %c is reserved for internal use as a gap, and cannot be used as a character in alphabet %s", c, alphabet.c_str());
    Assert (c != '>', "Character > is reserved by FASTA format as a delimiter, and cannot be used as a character in alphabet %s", alphabet.c_str());
    Assert (c != ' ' && c != '\t' && c != '\n', "Whitespace characters are reserved for internal use as delimiters, and cannot be used as alphabet characters");
    s.insert (c);
  }
  Assert (s.find(wild) == s.end(), "Wildcard character %c is also a character in alphabet %s", wild, alphabet.c_str());
}

UnvalidatedAlphTok AlphabetOwner::tokenize (char c) const {
  return ::tokenize (c, alphabet);
}

bool AlphabetOwner::isValidSymbol (char c) const {
  return tokenize(c) >= 0;
}

string AlphabetOwner::alphabetSymbol (AlphTok tok) const {
  return string (1, alphabet[tok]);
}

AlphTok AlphabetOwner::tokenizeOrDie (char c) const {
  UnvalidatedAlphTok tok = tokenize (c);
  if (tok >= 0 && tok < alphabetSize())
    return tok;
  Abort ("Character '%c' is not a member of alphabet '%s'", c, alphabet.c_str());
  return 0;
}

gsl_matrix* AlphabetOwner::newAlphabetMatrix() const {
  return gsl_matrix_alloc (alphabetSize(), alphabetSize());
}

gsl_vector* AlphabetOwner::newAlphabetVector() const {
  return gsl_vector_alloc (alphabetSize());
}

vguard<FastSeq> AlphabetOwner::convertWildcards (const vguard<FastSeq>& seqs) const {
  vguard<FastSeq> result (seqs);
  for (auto& seq: result)
    for (size_t n = 0; n < seq.seq.size(); ++n)
      if (seq.seq[n] == Alignment::wildcardChar)
	seq.seq[n] = wildcard;
  return result;
}

RateModel::RateModel()
{ }

RateModel::RateModel (const RateModel& model)
{
  *this = model;
}

RateModel::RateModel (const string& alph, int components, char wildcard) {
  initAlphabet (alph, wildcard);
  cptWeight = vguard<double> (components, 1. / (double) components);
  for (int c = 0; c < components; ++c) {
    auto ip = newAlphabetVector();
    auto sr = newAlphabetMatrix();
    insProb.push_back (ip);
    subRate.push_back (sr);
  }
}

RateModel::~RateModel()
{
  clear();
}

RateModel& RateModel::operator= (const RateModel& model) {
  clear();

  initAlphabet (model.alphabet, model.wildcard);
  copyIndelParams (model);
  cptWeight = model.cptWeight;
  for (int c = 0; c < model.components(); ++c) {
    auto ip = newAlphabetVector();
    auto sr = newAlphabetMatrix();
    CheckGsl (gsl_vector_memcpy (ip, model.insProb[c]));
    CheckGsl (gsl_matrix_memcpy (sr, model.subRate[c]));
    insProb.push_back (ip);
    subRate.push_back (sr);
  }
  return *this;
}

RateModel& RateModel::copyIndelParams (const RateModel& model) {
  insRate = model.insRate;
  delRate = model.delRate;
  insExtProb = model.insExtProb;
  delExtProb = model.delExtProb;
  return *this;
}

void RateModel::clear() {
  for (auto& sr: subRate)
    if (sr)
      gsl_matrix_free (sr);
  subRate.clear();

  for (auto& ip: insProb)
    if (ip)
      gsl_vector_free (ip);
  insProb.clear();

  cptWeight.clear();
}

void RateModel::init (const string& alph, char wild) {
  clear();
  initAlphabet (alph, wild);

  for (int c = 0; c < components(); ++c) {
    insProb.push_back (newAlphabetVector());
    subRate.push_back (newAlphabetMatrix());
    cptWeight.push_back (1. / components());
  }
}

void RateModel::read (const JsonValue& json) {
  Assert (subRate.empty(), "RateModel already initialized");

  readAlphabet (json);

  const JsonMap jm (json);
  insRate = jm.getNumber("insrate");
  insExtProb = jm.getNumber("insextprob");
  delRate = jm.getNumber("delrate");
  delExtProb = jm.getNumber("delextprob");

  if (jm.containsType("mixture",JSON_ARRAY)) {
    JsonValue jmix = jm.getType ("mixture", JSON_ARRAY);
    for (JsonIterator iter = begin(jmix); iter != end(jmix); ++iter) {
      const JsonMap jmcpt (iter->value);
      readComponent (jmcpt);
    }
  } else
    readComponent (jm);

  const double norm = accumulate (cptWeight.begin(), cptWeight.end(), 0.);
  for (auto& cw: cptWeight)
    cw /= norm;
}

void RateModel::readComponent (const JsonMap& jm) {
  gsl_matrix* sr = newAlphabetMatrix();
  gsl_matrix_set_zero (sr);

  const JsonMap rateMatrix = jm.getObject ("subrate");
  for (AlphTok i = 0; i < alphabet.size(); ++i) {
    const string si = alphabetSymbol(i);
    if (rateMatrix.contains (si)) {
      const JsonMap rateMatrixRow = rateMatrix.getObject (si);
      for (AlphTok j = 0; j < alphabet.size(); ++j) {
	if (j != i) {
	  const string sj = alphabetSymbol(j);
	  if (rateMatrixRow.contains (sj)) {
	    const double rate = rateMatrixRow.getNumber (sj);
	    *(gsl_matrix_ptr (sr, i, j)) += rate;
	    *(gsl_matrix_ptr (sr, i, i)) -= rate;
	  }
	}
      }
    }
  }

  gsl_vector* ip;
  if (jm.contains("rootprob")) {
    const JsonMap insVec = jm.getObject ("rootprob");
    ip = newAlphabetVector();
    gsl_vector_set_zero (ip);

    for (AlphTok i = 0; i < alphabet.size(); ++i) {
      const string si = alphabetSymbol(i);
      if (insVec.contains (si))
	gsl_vector_set (ip, i, insVec.getNumber(si));
    }
  } else
    ip = getEqmProbVector (sr);
  
  cptWeight.push_back (jm.containsType("weight",JSON_NUMBER) ? jm.getNumber("weight") : 1);
  insProb.push_back (ip);
  subRate.push_back (sr);
}

void RateModel::write (ostream& out) const {
  out << "{" << endl;
  out << " \"alphabet\": " << quoted_escaped(alphabet) << "," << endl;
  if (wildcard != Alignment::wildcardChar)
    out << " \"wildcard\": " << quoted_escaped(string(1,wildcard)) << "," << endl;
  out << " \"insrate\": " << insRate << "," << endl;
  out << " \"insextprob\": " << insExtProb << "," << endl;
  out << " \"delrate\": " << delRate << "," << endl;
  out << " \"delextprob\": " << delExtProb << "," << endl;
  if (components() > 1) {
    out << " \"mixture\": [" << endl;
    for (int c = 0; c < components(); ++c) {
      out << "  {" << endl;
      writeComponent (c, out);
      out << "  }" << (c < components() - 1 ? "," : "") << endl;
    }
    out << " ]" << endl;
  } else
    writeComponent (0, out);
  out << "}" << endl;
}

void RateModel::writeComponent (int cpt, ostream& out) const {
  const string indent (components() > 1 ? "   " : " ");
  if (components() > 1)
    out << indent << "\"weight\": " << cptWeight[cpt] << "," << endl;
  out << indent << "\"rootprob\":" << endl;
  out << indent << "{";
  for (AlphTok i = 0; i < alphabetSize(); ++i)
    out << (i>0 ? "," : "") << "\n" << indent << " \"" << alphabetSymbol(i) << "\": " << gsl_vector_get(insProb[cpt],i);
  out << endl;
  out << indent << "}," << endl;
  out << indent << "\"subrate\":" << endl;
  out << indent << "{" << endl;
  for (AlphTok i = 0; i < alphabetSize(); ++i) {
    out << indent << " \"" << alphabetSymbol(i) << "\": {";
    for (AlphTok j = 0; j < alphabetSize(); ++j)
      if (i != j)
	out << " \"" << alphabetSymbol(j) << "\": " << gsl_matrix_get(subRate[cpt],i,j) << (j < alphabetSize() - (i == alphabetSize() - 1 ? 2 : 1) ? "," : "");
    out << " }" << (i < alphabetSize() - 1 ? "," : "") << endl;
  }
  out << indent << "}" << endl;
}

gsl_vector* RateModel::getEqmProbVector (gsl_matrix* sr) {
  const size_t alphSize = sr->size1;
  // find eqm via QR decomposition
  gsl_matrix* QR = gsl_matrix_alloc (alphSize + 1, alphSize);
  for (AlphTok j = 0; j < alphSize; ++j) {
    for (AlphTok i = 0; i < alphSize; ++i)
      gsl_matrix_set (QR, i, j, gsl_matrix_get (sr, j, i));
    gsl_matrix_set (QR, alphSize, j, 1);
  }

  gsl_vector* tau = gsl_vector_alloc (alphSize);

  CheckGsl (gsl_linalg_QR_decomp (QR, tau));

  gsl_vector* b = gsl_vector_alloc (alphSize + 1);
  gsl_vector_set_zero (b);
  gsl_vector_set (b, alphSize, 1);

  gsl_vector* residual = gsl_vector_alloc (alphSize + 1);
  gsl_vector* eqm = gsl_vector_alloc (alphSize);
  
  CheckGsl (gsl_linalg_QR_lssolve (QR, tau, b, eqm, residual));

  // make sure no "probabilities" have become negative & (re-)normalize
  double eqmNorm = 0;
  for (AlphTok i = 0; i < alphSize; ++i) {
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

vguard<gsl_matrix*> RateModel::getSubProbMatrix (double t) const {
  vguard<gsl_matrix*> v;
  for (int c = 0; c < components(); ++c) {
    gsl_matrix* m = newAlphabetMatrix();
    gsl_matrix* rt = newAlphabetMatrix();
    CheckGsl (gsl_matrix_memcpy (rt, subRate[c]));
    CheckGsl (gsl_matrix_scale (rt, t));
    CheckGsl (gsl_linalg_exponential_ss (rt, m, GSL_PREC_DOUBLE));
    gsl_matrix_free (rt);
    v.push_back (m);
  }
  return v;
}

double RateModel::expectedSubstitutionRate() const {
  double R = 0;
  for (int c = 0; c < components(); ++c) {
    gsl_vector* eqm = getEqmProbVector (subRate[c]);
    for (AlphTok i = 0; i < alphabetSize(); ++i)
      for (AlphTok j = 0; j < alphabetSize(); ++j)
	if (i != j)
	  R += cptWeight[c] * gsl_vector_get (eqm, i) * gsl_matrix_get (subRate[c], i, j);
    gsl_vector_free (eqm);
  }
  return R;
}

RateModel RateModel::normalizeSubstitutionRate() const {
  return scaleRates (1. / expectedSubstitutionRate());
}

RateModel RateModel::scaleRates (double multiplier) const {
  return scaleRates (multiplier, multiplier);
}

RateModel RateModel::scaleRates (double substMultiplier, double indelMultiplier) const {
  RateModel model (*this);
  for (auto& sub: model.subRate)
    CheckGsl (gsl_matrix_scale (sub, substMultiplier));
  model.insRate *= indelMultiplier;
  model.delRate *= indelMultiplier;
  return model;
}

double RateModel::expectedInsertionLength() const {
  return 1. / (1. - insExtProb);
}

double RateModel::expectedDeletionLength() const {
  return 1. / (1. - delExtProb);
}

ProbModel::ProbModel (const RateModel& model, double t)
  : AlphabetOwner (model),
    t (t),
    ins (1 - exp (-model.insRate * t)),
    del (1 - exp (-model.delRate * t)),
    insExt (model.insExtProb),
    delExt (model.delExtProb),
    insWait (IndelCounts::decayWaitTime (model.insRate, t)),
    delWait (IndelCounts::decayWaitTime (model.delRate, t)),
    cptWeight (model.cptWeight),
    insVec (model.components()),
    subMat (model.getSubProbMatrix (t))
{
  for (int c = 0; c < model.components(); ++c) {
    insVec[c] = model.newAlphabetVector();
    CheckGsl (gsl_vector_memcpy (insVec[c], model.insProb[c]));
  }
}

ProbModel::~ProbModel() {
  for (auto& sr: subMat)
    gsl_matrix_free (sr);
  for (auto& iv: insVec)
    gsl_vector_free (iv);
}

double ProbModel::transProb (State src, State dest) const {
  switch (src) {
  case Match:
    switch (dest) {
    case Match:
      return (1 - ins) * (1 - del);
    case Insert:
      return ins;
    case Delete:
      return (1 - ins) * del;
    case End:
      return 1 - ins;
    default:
      break;
    }
    break;
  case Insert:
    switch (dest) {
    case Match:
      return (1 - insExt) * (1 - del);
    case Insert:
      return insExt;
    case Delete:
      return (1 - insExt) * del;
    case End:
      return 1 - insExt;
    default:
      break;
    }
    break;
  case Delete:
    switch (dest) {
    case Match:
    case End:
      return 1 - delExt;
    case Insert:
      return 0;
    case Delete:
      return delExt;
    default:
      break;
    }
    break;
  default:
    break;
  }
  return 0;
}

void ProbModel::write (ostream& out) const {
  out << "{" << endl;
  out << " \"alphabet\": \"" << alphabet << "\"," << endl;
  out << " \"insBegin\": " << ins << "," << endl;
  out << " \"insExtend\": " << insExt << "," << endl;
  out << " \"delBegin\": " << del << "," << endl;
  out << " \"delExtend\": " << delExt << "," << endl;

  if (components() > 1) {
    out << " \"mixture\": [" << endl;
    for (int c = 0; c < components(); ++c) {
      out << "  {" << endl;
      writeComponent (c, out);
      out << "  }" << (c < components() - 1 ? "," : "") << endl;
    }
    out << " ]" << endl;
  } else
    writeComponent (0, out);
  out << "}" << endl;
}

void ProbModel::writeComponent (int cpt, ostream& out) const {
  const string indent (components() > 1 ? "   " : " ");
  out << indent << "\"insVec\": {" << endl;
  for (AlphTok i = 0; i < alphabetSize(); ++i)
    out << indent << " \"" << alphabetSymbol(i) << "\": " << gsl_vector_get(insVec[cpt],i) << (i < alphabetSize() - 1 ? ",\n" : "");
  out << endl << indent << "}," << endl;
  out << indent << "\"subMat\": {" << endl;
  for (AlphTok i = 0; i < alphabetSize(); ++i) {
    out << indent << " \"" << alphabetSymbol(i) << "\": {" << endl;
    for (AlphTok j = 0; j < alphabetSize(); ++j)
      out << indent << "  \"" << alphabetSymbol(j) << "\": " << gsl_matrix_get(subMat[cpt],i,j) << (j < alphabetSize() - 1 ? ",\n" : "");
    out << endl << indent << " }" << (i < alphabetSize() - 1 ? "," : "") << endl;
  }
  out << indent << "}" << endl;
}

ProbModel::State ProbModel::getState (bool parentUngapped, bool childUngapped) {
  if (parentUngapped)
    return childUngapped ? Match : Delete;
  return childUngapped ? Insert : End;
}

LogProbModel::LogProbModel (const ProbModel& pm)
  : logCptWeight (pm.components()),
    logInsProb (pm.components(), vguard<LogProb> (pm.alphabetSize())),
    logSubProb (pm.components(), vguard<vguard<LogProb> > (pm.alphabetSize(), vguard<LogProb> (pm.alphabetSize())))
{
  for (int c = 0; c < components(); ++c) {
    logCptWeight[c] = log (pm.cptWeight[c]);
    logInsProb[c] = log_gsl_vector (pm.insVec[c]);
    for (AlphTok i = 0; i < pm.alphabetSize(); ++i)
      for (AlphTok j = 0; j < pm.alphabetSize(); ++j)
	logSubProb[c][i][j] = log (gsl_matrix_get (pm.subMat[c], i, j));
  }
}

double RateModel::mlDistance (const FastSeq& x, const FastSeq& y, int maxIterations) const {
  LogThisAt(7,"Estimating distance from " << x.name << " to " << y.name << endl);
  map<pair<AlphTok,AlphTok>,int> pairCount;
  Assert (x.length() == y.length(), "Sequences %s and %s have different lengths (%u, %u)", x.name.c_str(), y.name.c_str(), x.length(), y.length());
  for (size_t col = 0; col < x.length(); ++col) {
    const char ci = x.seq[col];
    const char cj = y.seq[col];
    if (!Alignment::isGap(ci) && !Alignment::isGap(cj) && !Alignment::isWildcard(ci) && !Alignment::isWildcard(cj)) {
      const UnvalidatedAlphTok toki = tokenize(ci), tokj = tokenize(cj);
      if (toki >= 0 && tokj >= 0)
	++pairCount[pair<AlphTok,AlphTok> (tokenizeOrDie(ci), tokenizeOrDie(cj))];
    }
  }
  if (LoggingThisAt(7)) {
    LogThisAt(7,"Counts:");
    for (const auto& pc : pairCount)
      LogThisAt(7," " << pc.second << "*" << alphabet[pc.first.first] << alphabet[pc.first.second]);
    LogThisAt(7,endl);
  }
  const DistanceMatrixParams dmp (pairCount, *this);
  const double t = dmp.tML (maxIterations);
  LogThisAt(6,"Distance from " << x.name << " to " << y.name << " is " << t << endl);
  return t;
}

vguard<vguard<double> > RateModel::distanceMatrix (const vguard<FastSeq>& gappedSeq, int maxIterations) const {
  vguard<vguard<double> > dist (gappedSeq.size(), vguard<double> (gappedSeq.size()));
  ProgressLog (plog, 4);
  plog.initProgress ("Distance matrix (%d rows)", gappedSeq.size());
  const size_t pairs = (gappedSeq.size() - 1) * gappedSeq.size() / 2;
  size_t n = 0;
  for (size_t i = 0; i + 1 < gappedSeq.size(); ++i)
    for (size_t j = i + 1; j < gappedSeq.size(); ++j) {
      plog.logProgress (n / (double) pairs, "computing entry %d/%d", n + 1, pairs);
      ++n;
      dist[i][j] = dist[j][i] = mlDistance (gappedSeq[i], gappedSeq[j], maxIterations);
    }
  if (LoggingThisAt(3)) {
    LogThisAt(3,"Distance matrix (" << dist.size() << " rows):" << endl);
    for (const auto& row : dist)
      LogThisAt(3,to_string_join(row) << endl);
  }
  return dist;
}

double distanceMatrixNegLogLike (double t, void *params) {
  const DistanceMatrixParams& dmp = *(const DistanceMatrixParams*) params;
  vguard<gsl_matrix*> sub = dmp.model.getSubProbMatrix (t);
  double ll = 0;
  for (auto& pc : dmp.pairCount) {
    double p = 0;
    for (int c = 0; c < dmp.model.components(); ++c)
      p += dmp.model.cptWeight[c] * gsl_matrix_get(sub[c],pc.first.first,pc.first.second);
    ll += log(p) * (double) pc.second;
  }
  for (auto& sr: sub)
    gsl_matrix_free (sr);
  return -ll;
}

double DistanceMatrixParams::negLogLike (double t) const {
  return distanceMatrixNegLogLike (t, (void*) this);
}

double DistanceMatrixParams::tJC() const {
  int same = 0, diff = 0;
  for (const auto& pc: pairCount)
    if (pc.first.first == pc.first.second)
      same += pc.second;
    else
      diff += pc.second;
  const double pDiff = diff / (double) (same + diff);
  const double A = (double) model.alphabetSize();
  if (pDiff >= (A - 1) / A)
    return numeric_limits<double>::infinity();
  return -((A-1) / A) * log (1 - (A/(A-1)) * pDiff) / model.expectedSubstitutionRate();
}

double DistanceMatrixParams::tML (int maxIterations) const {
  const double tMin = 1e-9, tMax = 10;
  const double tjc = min (tMax, max (tMin, tJC()));
  if (maxIterations <= 0)
    return tjc;
  
  const gsl_min_fminimizer_type *T;
  gsl_min_fminimizer *s;

  gsl_function F;
  F.function = &distanceMatrixNegLogLike;
  F.params = (void*) this;

  T = gsl_min_fminimizer_goldensection; // gsl_min_fminimizer_brent;
  s = gsl_min_fminimizer_alloc (T);

  double t;
  const double tLower = min (tMin, tjc/2), tUpper = max (tMax, tjc*2);
  const double llLower = negLogLike(tLower), llUpper = negLogLike(tUpper);
  const double lljc = negLogLike(tjc);
  LogThisAt(8,"tML: f(" << tLower << ") = " << llLower << ", f(" << tjc << ") = " << lljc << ", f(" << tUpper << ") = " << llUpper << endl);

  if (lljc < llLower && lljc < llUpper)
    t = tjc;
  else {
    bool foundGuess = false;
    double tScanLower = tLower, tScanUpper = tUpper;
    const double nScanSteps = 4;
    while (!foundGuess && tScanUpper - tScanLower > tLower) {
      const double step = (tScanUpper - tScanLower) / nScanSteps;
      LogThisAt(9,"tML: Scanning from " << tScanLower << " to " << tScanUpper << " step " << step << endl);
      for (double x = tScanLower; x < tScanUpper && !foundGuess; x += step) {
	const double ll = negLogLike(x);
	LogThisAt(8,"tML: f(" << x << ") = " << ll << endl);
	if (ll < llLower && ll < llUpper) {
	  foundGuess = true;
	  t = x;
	}
      }
      if (!foundGuess) {
	if (llLower < llUpper)
	  tScanUpper = (tScanLower + tScanUpper) / 2;
	else
	  tScanLower = (tScanLower + tScanUpper) / 2;
      }
    }
    if (!foundGuess)
      return llLower < llUpper ? tLower : tUpper;
  }

  LogThisAt(8,"Initializing with t=" << t << ", tLower=" << tLower << ", tUpper=" << tUpper << endl);
  gsl_min_fminimizer_set (s, &F, t, tLower, tUpper);

  const double tML_convergence = .01; // converge to 1%
  for (int iter = 0; iter < maxIterations; ++iter) {
    (void) gsl_min_fminimizer_iterate (s);

    t = gsl_min_fminimizer_x_minimum (s);
    const double a = gsl_min_fminimizer_x_lower (s);
    const double b = gsl_min_fminimizer_x_upper (s);

    const int status = gsl_min_test_interval (a, b, 0., tML_convergence);
    LogThisAt(7,"Iteration #" << iter+1 << " tML: f(" << t << ") = " << negLogLike(t) << " (a=" << a << ", b=" << b << ", status=" << status << ")" << endl);
    
    if (status == GSL_SUCCESS)
      break;
  }

  gsl_min_fminimizer_free (s);

  return t;
}

void AlphabetOwner::writeSubCounts (ostream& out, const vguard<vguard<double> >& rootCounts, const vguard<vguard<vguard<double> > >& subCountsAndWaitTimes, size_t indent) const {
  const string ind (indent, ' ');
  if (rootCounts.size() > 1) {
    out << ind << "{" << endl;
    out << ind << " \"mixture\": [" << endl;
    for (int c = 0; c < (int) rootCounts.size(); ++c) {
      writeSubCountsComponent (out, rootCounts[c], subCountsAndWaitTimes[c], indent + 2);
      if (c < rootCounts.size() - 1)
	out << ",";
      out << endl;
    }
    out << ind << " ]" << endl;
    out << ind << "}";
  } else
    writeSubCountsComponent (out, rootCounts[0], subCountsAndWaitTimes[0], indent);
}

void AlphabetOwner::writeSubCountsComponent (ostream& out, const vguard<double>& rootCounts, const vguard<vguard<double> >& subCountsAndWaitTimes, size_t indent) const {
  const string ind (indent, ' ');
  out << ind << "{" << endl;
  out << ind << " \"root\":" << endl;
  out << ind << "  {";
  for (AlphTok i = 0; i < alphabetSize(); ++i)
    out << (i == 0 ? "" : ",") << endl << ind << "   \"" << alphabetSymbol(i) << "\": " << rootCounts[i];
  out << endl << ind << "  }," << endl;
  out << ind << " \"sub\":" << endl;
  out << ind << "  {";
  for (AlphTok i = 0; i < alphabetSize(); ++i) {
    out << (i == 0 ? "" : ",") << endl << ind << "   \"" << alphabetSymbol(i) << "\": {";
    for (AlphTok j = 0; j < alphabetSize(); ++j)
      if (i != j)
	out << (j == (i == 0 ? 1 : 0) ? "" : ",") << " " << "\"" << alphabetSymbol(j) << "\": " << subCountsAndWaitTimes[i][j];
    out << " }";
  }
  out << endl;
  out << ind << "  }," << endl;
  out << ind << " \"wait\":" << endl;
  out << ind << "  {";
  for (AlphTok i = 0; i < alphabetSize(); ++i)
    out << (i == 0 ? "" : ",") << endl << ind << "   \"" << alphabetSymbol(i) << "\": " << subCountsAndWaitTimes[i][i];
  out << endl;
  out << ind << "  }" << endl;
  out << ind << "}";
}

IndelCounts::IndelCounts (double pseudocount, double pseudotime)
  : ins(pseudocount),
    del(pseudocount),
    insExt(pseudocount),
    delExt(pseudocount),
    insTime(pseudotime),
    delTime(pseudotime),
    lp(0)
{ }

IndelCounts& IndelCounts::operator+= (const IndelCounts& c) {
  ins += c.ins;
  del += c.del;
  insExt += c.insExt;
  delExt += c.delExt;
  insTime += c.insTime;
  delTime += c.delTime;
  lp += c.lp;
  return *this;
}

IndelCounts& IndelCounts::operator*= (double w) {
  ins *= w;
  del *= w;
  insExt *= w;
  delExt *= w;
  insTime *= w;
  delTime *= w;
  lp *= w;
  return *this;
}

IndelCounts IndelCounts::operator+ (const IndelCounts& c) const {
  IndelCounts result (*this);
  result += c;
  return result;
}

IndelCounts IndelCounts::operator* (double w) const {
  IndelCounts result (*this);
  result *= w;
  return result;
}

EigenCounts::EigenCounts()
{ }

EigenCounts::EigenCounts (size_t components, size_t alphabetSize)
  : indelCounts(),
    rootCount (components, vguard<double> (alphabetSize, 0)),
    eigenCount (components, vguard<vguard<gsl_complex> > (alphabetSize, vguard<gsl_complex> (alphabetSize, gsl_complex_rect(0,0))))
{ }

EigenCounts& EigenCounts::operator+= (const EigenCounts& c) {
  indelCounts += c.indelCounts;

  if (components() == 0) {
    rootCount = c.rootCount;
    eigenCount = c.eigenCount;
  } else if (c.components() > 0) {
    for (int cpt = 0; cpt < components(); ++cpt) {
      for (size_t n = 0; n < rootCount[cpt].size(); ++n)
	rootCount[cpt][n] += c.rootCount.at(cpt).at(n);

      for (size_t i = 0; i < eigenCount[cpt].size(); ++i)
	for (size_t j = 0; j < eigenCount[cpt][i].size(); ++j)
	  eigenCount[cpt][i][j] = gsl_complex_add (eigenCount[cpt][i][j], c.eigenCount.at(cpt).at(i).at(j));
    }
  }

  return *this;
}

EigenCounts& EigenCounts::operator*= (double w) {
  indelCounts *= w;
  for (int cpt = 0; cpt < components(); ++cpt) {
    for (auto& c : rootCount[cpt])
      c *= w;
    for (auto& sc : eigenCount[cpt])
      for (auto& c : sc)
	c = gsl_complex_mul_real (c, w);
  }
  return *this;
}

EigenCounts EigenCounts::operator+ (const EigenCounts& c) const {
  EigenCounts result (*this);
  result += c;
  return result;
}

EigenCounts EigenCounts::operator* (double w) const {
  EigenCounts result (*this);
  result *= w;
  return result;
}

EventCounts::EventCounts (const AlphabetOwner& alph, int components, double pseudo)
  : AlphabetOwner (alph),
    indelCounts (pseudo, pseudo),
    rootCount (components, vguard<double> (alph.alphabetSize(), pseudo)),
    subCount (components, vguard<vguard<double> > (alph.alphabetSize(), vguard<double> (alph.alphabetSize(), pseudo)))
{ }

EventCounts& EventCounts::operator+= (const EventCounts& c) {
  Assert (alphabet == c.alphabet, "Alphabets don't match");

  indelCounts += c.indelCounts;

  for (int cpt = 0; cpt < components(); ++cpt) {
    for (size_t n = 0; n < rootCount[cpt].size(); ++n)
      rootCount[cpt][n] += c.rootCount.at(cpt).at(n);

    for (size_t i = 0; i < subCount[cpt].size(); ++i)
      for (size_t j = 0; j < subCount[cpt][i].size(); ++j)
	subCount[cpt][i][j] += c.subCount.at(cpt).at(i).at(j);
  }

  return *this;
}

EventCounts& EventCounts::operator*= (double w) {
  indelCounts *= w;
  for (int cpt = 0; cpt < components(); ++cpt) {
    for (auto& c : rootCount[cpt])
      c *= w;
    for (auto& sc : subCount[cpt])
      for (auto& c : sc)
	c *= w;
  }
  return *this;
}

EventCounts EventCounts::operator+ (const EventCounts& c) const {
  EventCounts result (*this);
  result += c;
  return result;
}

EventCounts EventCounts::operator* (double w) const {
  EventCounts result (*this);
  result *= w;
  return result;
}

void IndelCounts::accumulateIndelCounts (const RateModel& model, double time, const AlignRowPath& parent, const AlignRowPath& child, double weight) {
  const double insWait = decayWaitTime (model.insRate, time), delWait = decayWaitTime (model.delRate, time);
  ProbModel pm (model, time);
  ProbModel::State state, next;
  state = ProbModel::Match;
  for (size_t col = 0; col < parent.size(); ++col) {
    if (parent[col] && child[col])
      next = ProbModel::Match;
    else if (parent[col])
      next = ProbModel::Delete;
    else if (child[col])
      next = ProbModel::Insert;
    else
      continue;
    switch (next) {
    case ProbModel::Match:
      if (state == next) {
	insTime += weight * time;
	delTime += weight * time;
      }
      break;
    case ProbModel::Insert:
      if (state == next)
	insExt += weight;
      else {
	ins += weight;
	insTime += weight * insWait;
      }
      break;
    case ProbModel::Delete:
      if (state == next)
	delExt += weight;
      else {
	del += weight;
	delTime += weight * delWait;
      }
      break;
    default:
      break;
    }
    lp += log (pm.transProb (state, next)) * weight;
    state = next;
  }
  lp += log (pm.transProb (state, ProbModel::End)) * weight;

  LogThisAt(9,"Accumulated indel counts on branch of length " << time << " (weight=" << weight << "):\n" << toJson());
}

void IndelCounts::accumulateIndelCounts (const RateModel& model, const Tree& tree, const AlignPath& align, double weight) {
  for (TreeNodeIndex node = 0; node < tree.nodes() - 1; ++node)
    accumulateIndelCounts (model, tree.branchLength(node), align.at(tree.parentNode(node)), align.at(node), weight);
}

void EigenCounts::accumulateSubstitutionCounts (const RateModel& model, const Tree& tree, const vguard<FastSeq>& gapped, double weight) {
  AlignColSumProduct colSumProd (model, tree, gapped);

  EigenCounts c (model.components(), model.alphabetSize());

  while (!colSumProd.alignmentDone()) {
    colSumProd.fillUp();
    colSumProd.fillDown();
    colSumProd.accumulateEigenCounts (c.rootCount, c.eigenCount);
    c.indelCounts.lp += colSumProd.columnLogLikelihood();
    colSumProd.nextColumn();
  }

  c *= weight;
  *this += c;
}

void EigenCounts::accumulateCounts (const RateModel& model, const Alignment& align, const Tree& tree, bool updateIndelCounts, bool updateSubstCounts, double weight) {
  if (updateIndelCounts)
    indelCounts.accumulateIndelCounts (model, tree, align.path, weight);
  if (updateSubstCounts)
    accumulateSubstitutionCounts (model, tree, align.gapped(), weight);
}

EventCounts EigenCounts::transform (const RateModel& model) const {
  EventCounts c (model, model.components());
  EigenModel eigen (model);
  c.indelCounts = indelCounts;
  c.rootCount = rootCount;
  c.subCount = eigen.getSubCounts (eigenCount);
  return c;
}

void EventCounts::writeJson (ostream& out) const {
  out << "{" << endl;
  out << " \"alphabet\": \"" << alphabet << "\"," << endl;
  out << " \"indel\":" << endl;
  indelCounts.writeJson (out, 2);
  out << "," << endl;
  out << " \"sub\":" << endl;
  writeSubCounts (out, rootCount, subCount, 2);
  out << "," << endl;
  out << " \"logLikelihood\": " << indelCounts.lp << endl;
  out << "}" << endl;
}

void IndelCounts::writeJson (ostream& out, const size_t indent) const {
  const string ind (indent, ' ');
  out << ind << "{" << endl;
  out << ind << " \"ins\": " << ins << "," << endl;
  out << ind << " \"del\": " << del << "," << endl;
  out << ind << " \"insExt\": " << insExt << "," << endl;
  out << ind << " \"delExt\": " << delExt << "," << endl;
  out << ind << " \"insTime\": " << insTime << endl;
  out << ind << " \"delTime\": " << delTime << endl;
  out << ind << "}";
}

void IndelCounts::read (const JsonValue& json) {
  JsonMap jm (json);
  ins = jm.getNumber("ins");
  del = jm.getNumber("del");
  insExt = jm.getNumber("insExt");
  delExt = jm.getNumber("delExt");
  insTime = jm.getNumber("insTime");
  delTime = jm.getNumber("delTime");
}

string IndelCounts::toJson() const {
  ostringstream out;
  writeJson (out);
  out << endl;
  return out.str();
}

string EventCounts::toJson() const {
  ostringstream out;
  writeJson (out);
  out << endl;
  return out.str();
}

void EventCounts::read (const JsonValue& json) {
  readAlphabet (json);
  rootCount.clear();
  subCount.clear();
  JsonMap jm (json);
  indelCounts.read (jm.getType("indel",JSON_OBJECT));
  indelCounts.lp = jm.getNumber("logLikelihood");

  if (jm.containsType("mixture",JSON_ARRAY)) {
    JsonValue jmix = jm.getType ("mixture", JSON_ARRAY);
    for (JsonIterator iter = begin(jmix); iter != end(jmix); ++iter) {
      const JsonMap jmcpt (iter->value);
      readComponent (jmcpt);
    }
  } else
    readComponent (jm);
}

void EventCounts::readComponent (const JsonMap& jm) {
  vguard<double> rc (alphabetSize(), 0.);
  vguard<vguard<double> > sc (alphabetSize(), vguard<double> (alphabetSize(), 0.));
  JsonMap subBlock = jm.getObject("sub");
  JsonMap root = subBlock.getObject("root");
  JsonMap sub = subBlock.getObject("sub");
  JsonMap wait = subBlock.getObject("wait");
  for (AlphTok i = 0; i < alphabetSize(); ++i) {
    const string si = alphabetSymbol(i);
    rc[i] = root.getNumber(si);
    sc[i][i] = wait.getNumber(si);
    JsonMap sub_i = sub.getObject(si);
    for (AlphTok j = 0; j < alphabetSize(); ++j)
      if (i != j) {
	const string sj = alphabetSymbol(j);
	sc[i][j] = sub_i.getNumber(sj);
      }
  }
  rootCount.push_back (rc);
  subCount.push_back (sc);
}

void EventCounts::optimize (RateModel& model, bool fitIndelRates, bool fitSubstRates) const {
  if (model.alphabet != alphabet)
    model.init (alphabet);
  
  LogThisAt(9,"Optimizing model for the following expected counts:\n" << toJson());
  
  if (fitSubstRates) {
    vguard<double> cptCount;
    double totalCptCount = 0;
    for (int cpt = 0; cpt < components(); ++cpt) {
      const double insNorm = accumulate (rootCount[cpt].begin(), rootCount[cpt].end(), 0.);
      for (AlphTok i = 0; i < alphabetSize(); ++i)
	gsl_vector_set (model.insProb[cpt], i, rootCount[cpt][i] / insNorm);

      for (AlphTok i = 0; i < alphabetSize(); ++i) {
	double r_ii = 0;
	for (AlphTok j = 0; j < alphabetSize(); ++j)
	  if (j != i) {
	    const double r_ij = subCount[cpt][i][j] / subCount[cpt][i][i];
	    gsl_matrix_set (model.subRate[cpt], i, j, r_ij);
	    r_ii -= r_ij;
	  }
	gsl_matrix_set (model.subRate[cpt], i, i, r_ii);
      }
      cptCount.push_back (insNorm);
      totalCptCount += insNorm;
    }
    for (int cpt = 0; cpt < components(); ++cpt)
      model.cptWeight[cpt] = cptCount[cpt] / totalCptCount;
  }

  if (fitIndelRates) {
    model.insRate = indelCounts.ins / indelCounts.insTime;
    model.delRate = indelCounts.del / indelCounts.delTime;
    model.insExtProb = indelCounts.insExt / (indelCounts.insExt + indelCounts.ins);
    model.delExtProb = indelCounts.delExt / (indelCounts.delExt + indelCounts.del);
  }
}

double EventCounts::logPrior (const RateModel& model, bool includeIndelRates, bool includeSubstRates) const {
  double lp = 0;
  if (includeIndelRates)
    lp += logGammaPdf (model.insRate, indelCounts.ins, indelCounts.insTime)
      + logGammaPdf (model.delRate, indelCounts.del, indelCounts.delTime)
      + logBetaPdf (model.insExtProb, indelCounts.insExt, indelCounts.ins)
      + logBetaPdf (model.delExtProb, indelCounts.delExt, indelCounts.del);
  if (includeSubstRates)
    for (int cpt = 0; cpt < components(); ++cpt) {
      lp += logDirichletPdf (gsl_vector_to_stl(model.insProb[cpt]), rootCount[cpt]);
      for (AlphTok i = 0; i < alphabetSize(); ++i)
	for (AlphTok j = 0; j < alphabetSize(); ++j)
	  if (i != j)
	    lp += logGammaPdf (gsl_matrix_get (model.subRate[cpt], i, j), subCount[cpt][i][j], subCount[cpt][i][i]);
    }
  return lp;
}

double xlogy (double x, double y) {
  return x > 0 && y > 0 ? x*log(y) : 0;
}

double EventCounts::expectedLogLikelihood (const RateModel& model) const {
  double lp = -model.insRate * indelCounts.insTime
    + xlogy (indelCounts.ins, model.insRate)
    -model.delRate * indelCounts.delTime
    + xlogy (indelCounts.del, model.delRate)
    + xlogy (indelCounts.insExt, model.insExtProb)
    + xlogy (indelCounts.ins, 1 - model.insExtProb)
    + xlogy (indelCounts.delExt, model.delExtProb)
    + xlogy (indelCounts.del, 1 - model.delExtProb);
  
  for (int cpt = 0; cpt < components(); ++cpt)
    for (AlphTok i = 0; i < alphabetSize(); ++i) {
      const double exit_i = -gsl_matrix_get (model.subRate[cpt], i, i);
      lp += xlogy (rootCount[cpt][i], gsl_vector_get (model.insProb[cpt], i));
      lp -= exit_i * subCount[cpt][i][i];
      for (AlphTok j = 0; j < alphabetSize(); ++j)
	if (i != j)
	  lp += xlogy (subCount[cpt][i][j], gsl_matrix_get (model.subRate[cpt], i, j));
    }

  return lp;
}

double IndelCounts::decayWaitTime (double decayRate, double timeInterval) {
  return 1 / decayRate - timeInterval / (exp (decayRate*timeInterval) - 1);
}

EigenModel::EigenModel (const EigenModel& eigen)
  : model (eigen.model),
    eval (eigen.components()),
    evec (eigen.components()),
    evecInv (eigen.components()),
    ev (eigen.ev),
    ev_t (eigen.ev_t),
    exp_ev_t (eigen.exp_ev_t),
    isReal (eigen.isReal),
    realEval (eigen.realEval),
    realEvec (eigen.realEvec),
    realEvecInv (eigen.realEvecInv),
    real_ev_t (eigen.real_ev_t),
    real_exp_ev_t (eigen.real_exp_ev_t)
{
  for (int cpt = 0; cpt < eigen.components(); ++cpt) {
    eval[cpt] = gsl_vector_complex_alloc (eigen.model.alphabetSize());
    evec[cpt] = gsl_matrix_complex_alloc (eigen.model.alphabetSize(), eigen.model.alphabetSize());
    evecInv[cpt] = gsl_matrix_complex_alloc (eigen.model.alphabetSize(), eigen.model.alphabetSize());
    CheckGsl (gsl_vector_complex_memcpy (eval[cpt], eigen.eval[cpt]));
    CheckGsl (gsl_matrix_complex_memcpy (evec[cpt], eigen.evec[cpt]));
    CheckGsl (gsl_matrix_complex_memcpy (evecInv[cpt], eigen.evecInv[cpt]));
  }
}

EigenModel::EigenModel (const RateModel& model)
  : model (model),
    eval (model.components()),
    evec (model.components()),
    evecInv (model.components()),
    ev (model.components(), vguard<gsl_complex> (model.alphabetSize())),
    ev_t (model.components(), vguard<gsl_complex> (model.alphabetSize())),
    exp_ev_t (model.components(), vguard<gsl_complex> (model.alphabetSize())),
    isReal (model.components(), false),
    realEval (model.components(), vguard<double> (model.alphabetSize())),
    realEvec (model.components(), vguard<vguard<double> > (model.alphabetSize(), vguard<double> (model.alphabetSize()))),
    realEvecInv (model.components(), vguard<vguard<double> > (model.alphabetSize(), vguard<double> (model.alphabetSize()))),
    real_ev_t (model.components(), vguard<double> (model.alphabetSize())),
    real_exp_ev_t (model.components(), vguard<double> (model.alphabetSize()))
{
  for (int cpt = 0; cpt < model.components(); ++cpt) {
    eval[cpt] = gsl_vector_complex_alloc (model.alphabetSize());
    evec[cpt] = gsl_matrix_complex_alloc (model.alphabetSize(), model.alphabetSize());
    evecInv[cpt] = gsl_matrix_complex_alloc (model.alphabetSize(), model.alphabetSize());

    gsl_matrix *R = gsl_matrix_alloc (model.alphabetSize(), model.alphabetSize());
    gsl_matrix_memcpy (R, model.subRate[cpt]);
  
    gsl_eigen_nonsymmv_workspace *workspace = gsl_eigen_nonsymmv_alloc (model.alphabetSize());
    CheckGsl (gsl_eigen_nonsymmv (R, eval[cpt], evec[cpt], workspace));
    gsl_eigen_nonsymmv_free (workspace);
    gsl_matrix_free (R);

    gsl_matrix_complex *LU = gsl_matrix_complex_alloc (model.alphabetSize(), model.alphabetSize());
    gsl_permutation *perm = gsl_permutation_alloc (model.alphabetSize());
    int permSig = 0;
    gsl_matrix_complex_memcpy (LU, evec[cpt]);
    CheckGsl (gsl_linalg_complex_LU_decomp (LU, perm, &permSig));
    CheckGsl (gsl_linalg_complex_LU_invert (LU, perm, evecInv[cpt]));
    gsl_matrix_complex_free (LU);
    gsl_permutation_free (perm);

    for (AlphTok i = 0; i < model.alphabetSize(); ++i)
      ev[cpt][i] = gsl_vector_complex_get (eval[cpt], i);

    isReal[cpt] = true;
    for (AlphTok i = 0; isReal[cpt] && i < model.alphabetSize(); ++i) {
      isReal[cpt] = isReal[cpt] && EIGENMODEL_NEAR_REAL(ev[cpt][i]);
      for (AlphTok j = 0; isReal[cpt] && j < model.alphabetSize(); ++j)
	isReal[cpt] = isReal[cpt] && EIGENMODEL_NEAR_REAL(gsl_matrix_complex_get(evec[cpt],i,j))
	  && EIGENMODEL_NEAR_REAL(gsl_matrix_complex_get(evecInv[cpt],i,j));
    }

    if (isReal[cpt])
      for (AlphTok i = 0; isReal[cpt] && i < model.alphabetSize(); ++i) {
	realEval[cpt][i] = GSL_REAL (ev[cpt][i]);
	for (AlphTok j = 0; isReal[cpt] && j < model.alphabetSize(); ++j) {
	  realEvec[cpt][i][j] = GSL_REAL (gsl_matrix_complex_get (evec[cpt], i, j));
	  realEvecInv[cpt][i][j] = GSL_REAL (gsl_matrix_complex_get (evecInv[cpt], i, j));
	}
      }
  
    LogThisAt(8,"Component #" << cpt << endl
	      << "Eigenvalues:" << complexVectorToString(ev[cpt]) << endl
	      << "Right eigenvector matrix, V:" << endl << complexMatrixToString(evec[cpt])
	      << "Left eigenvector matrix, V^{-1}:" << endl << complexMatrixToString(evecInv[cpt])
	      << "Product V^{-1} * V:" << endl << tempComplexMatrixToString (evecInv_evec(cpt))
	      << "Reconstituted rate matrix:" << endl << tempComplexMatrixToString (getRateMatrix(cpt)));
  }
}

EigenModel::~EigenModel() {
  for (auto& e: eval)
    if (e)
      gsl_vector_complex_free (e);
  for (auto& e: evec)
    if (e)
      gsl_matrix_complex_free (e);
  for (auto& e: evecInv)
    if (e)
      gsl_matrix_complex_free (e);
}

void EigenModel::compute_exp_ev_t (double t, bool forceComplex) {
  for (int cpt = 0; cpt < model.components(); ++cpt) {
    if (isReal[cpt] && !forceComplex) {
      for (AlphTok i = 0; i < model.alphabetSize(); ++i) {
	real_ev_t[cpt][i] = realEval[cpt][i] * t;
	real_exp_ev_t[cpt][i] = exp (real_ev_t[cpt][i]);
      }
      LogThisAt(9,"Component #" << cpt << " exp(eigenvalue*" << t << "):" << join(real_exp_ev_t[cpt]," ") << endl);
    } else {
      for (AlphTok i = 0; i < model.alphabetSize(); ++i) {
	ev_t[cpt][i] = gsl_complex_mul_real (ev[cpt][i], t);
	exp_ev_t[cpt][i] = gsl_complex_exp (ev_t[cpt][i]);
      }
      LogThisAt(9,"Component #" << cpt << " exp(eigenvalue*" << t << "):" << complexVectorToString(exp_ev_t[cpt]));
    }
  }
}

gsl_matrix_complex* EigenModel::getRateMatrix (int cpt) const {
  gsl_matrix_complex* r = gsl_matrix_complex_alloc (model.alphabetSize(), model.alphabetSize());
  for (AlphTok i = 0; i < model.alphabetSize(); ++i)
    for (AlphTok j = 0; j < model.alphabetSize(); ++j) {
      gsl_complex rij = gsl_complex_rect (0, 0);
      for (AlphTok k = 0; k < model.alphabetSize(); ++k)
	rij = gsl_complex_add
	  (rij,
	   gsl_complex_mul (gsl_complex_mul (gsl_matrix_complex_get (evec[cpt], i, k),
					     gsl_matrix_complex_get (evecInv[cpt], k, j)),
			    ev[cpt][k]));
      gsl_matrix_complex_set (r, i, j, rij);
    }
  return r;
}

gsl_matrix_complex* EigenModel::evecInv_evec (int cpt) const {
  gsl_matrix_complex* e = gsl_matrix_complex_alloc (model.alphabetSize(), model.alphabetSize());
  for (AlphTok i = 0; i < model.alphabetSize(); ++i)
    for (AlphTok j = 0; j < model.alphabetSize(); ++j) {
      gsl_complex eij = gsl_complex_rect (0, 0);
      for (AlphTok k = 0; k < model.alphabetSize(); ++k)
	eij = gsl_complex_add
	  (eij,
	   gsl_complex_mul (gsl_matrix_complex_get (evec[cpt], i, k),
			    gsl_matrix_complex_get (evecInv[cpt], k, j)));
      gsl_matrix_complex_set (e, i, j, eij);
    }
  return e;
}

double EigenModel::getSubProb (int cpt, double t, AlphTok i, AlphTok j) const {
  ((EigenModel&) *this).compute_exp_ev_t (t);
  return getSubProbInner (cpt, t, i, j);
}

double EigenModel::getSubProbInner (int cpt, double t, AlphTok i, AlphTok j) const {
  if (isReal[cpt]) {
    double p = 0;
    for (AlphTok k = 0; k < model.alphabetSize(); ++k)
      p += realEvec[cpt][i][k] * realEvecInv[cpt][k][j] * real_exp_ev_t[cpt][k];
    return min (1., max (0., p));
  }
  gsl_complex p = gsl_complex_rect (0, 0);
  for (AlphTok k = 0; k < model.alphabetSize(); ++k)
    p = gsl_complex_add
      (p,
       gsl_complex_mul (gsl_complex_mul (gsl_matrix_complex_get (evec[cpt], i, k),
					 gsl_matrix_complex_get (evecInv[cpt], k, j)),
			exp_ev_t[cpt][k]));
  Assert (EIGENMODEL_NEAR_REAL(p), "Probability has imaginary part: p=(%g,%g)", GSL_REAL(p), GSL_IMAG(p));
  return min (1., max (0., GSL_REAL(p)));
}

vguard<gsl_matrix*> EigenModel::getSubProbMatrix (double t) const {
  vguard<gsl_matrix*> v;
  for (int cpt = 0; cpt < model.components(); ++cpt) {
    gsl_matrix* sub = gsl_matrix_alloc (model.alphabetSize(), model.alphabetSize());
    ((EigenModel&) *this).compute_exp_ev_t (t);
    for (AlphTok i = 0; i < model.alphabetSize(); ++i)
      for (AlphTok j = 0; j < model.alphabetSize(); ++j)
	gsl_matrix_set (sub, i, j, getSubProbInner (cpt, t, i, j));
    v.push_back (sub);
  }
  return v;
}

double EigenModel::getSubCount (int cpt, AlphTok a, AlphTok b, AlphTok i, AlphTok j, const gsl_matrix* sub, const gsl_matrix_complex* eSubCount) const {
  const double p_ab = gsl_matrix_get (sub, a, b);
  const double r_ij = gsl_matrix_get (model.subRate[cpt], i, j);
  gsl_complex c_ij = gsl_complex_rect (0, 0);
  for (AlphTok k = 0; k < model.alphabetSize(); ++k) {
    gsl_complex c_ijk = gsl_complex_rect (0, 0);
    for (AlphTok l = 0; l < model.alphabetSize(); ++l)
      c_ijk = gsl_complex_add
	(c_ijk,
	 gsl_complex_mul
	 (gsl_complex_mul
	  (gsl_matrix_complex_get (evec[cpt], j, l),
	   gsl_matrix_complex_get (evecInv[cpt], l, b)),
	  gsl_matrix_complex_get (eSubCount, k, l)));
    c_ij = gsl_complex_add (c_ij,
			    gsl_complex_mul
			    (gsl_complex_mul
			     (gsl_matrix_complex_get (evec[cpt], a, k),
			      gsl_matrix_complex_get (evecInv[cpt], k, i)),
			     c_ijk));
  }
  Assert (EIGENMODEL_NEAR_REAL(c_ij), "Count has imaginary part: c=(%g,%g)", GSL_REAL(c_ij), GSL_IMAG(c_ij));
  return max (0., (i == j ? 1. : r_ij) * GSL_REAL(c_ij) / p_ab);
}

void EigenModel::accumSubCounts (int cpt, vguard<vguard<double> >& count, AlphTok a, AlphTok b, double weight, const gsl_matrix* sub, const gsl_matrix_complex* eSubCount) const {
  for (AlphTok i = 0; i < model.alphabetSize(); ++i)
    for (AlphTok j = 0; j < model.alphabetSize(); ++j)
      count[i][j] += getSubCount (cpt, a, b, i, j, sub, eSubCount) * weight;
}

vguard<gsl_matrix_complex*> EigenModel::eigenSubCount (double t) const {
  vguard<gsl_matrix_complex*> v;
  ((EigenModel&) *this).compute_exp_ev_t (t, true);
  for (int cpt = 0; cpt < components(); ++cpt) {
    gsl_matrix_complex* esub = gsl_matrix_complex_alloc (model.alphabetSize(), model.alphabetSize());
    for (AlphTok i = 0; i < model.alphabetSize(); ++i)
      for (AlphTok j = 0; j < model.alphabetSize(); ++j) {
	const bool ev_eq = i == j || EIGENMODEL_NEAR_EQ_COMPLEX (ev[cpt][i], ev[cpt][j]);
	gsl_matrix_complex_set
	  (esub, i, j,
	   ev_eq
	   ? gsl_complex_mul_real (exp_ev_t[cpt][i], t)
	   : gsl_complex_div (gsl_complex_sub (exp_ev_t[cpt][i], exp_ev_t[cpt][j]),
			      gsl_complex_sub (ev[cpt][i], ev[cpt][j])));
      }

    LogThisAt(8,endl << "Component #" << cpt << " eigensubstitution matrix at time t=" << t << ":" << endl << complexMatrixToString(esub));
    v.push_back (esub);
  }

  return v;
}

vguard<vguard<vguard<double> > > EigenModel::getSubCounts (const vguard<vguard<vguard<gsl_complex> > >& eigenCounts) const {
  vguard<vguard<vguard<double> > > v;
  for (int cpt = 0; cpt < components(); ++cpt) {
    LogThisAt(8,"Component #" << cpt << " eigencounts matrix:" << endl << complexMatrixToString(eigenCounts[cpt]) << endl);
    vguard<vguard<double> > counts (model.alphabetSize(), vguard<double> (model.alphabetSize(), 0));
    for (AlphTok i = 0; i < model.alphabetSize(); ++i)
      for (AlphTok j = 0; j < model.alphabetSize(); ++j) {
	gsl_complex c = gsl_complex_rect (0, 0);
	for (AlphTok k = 0; k < model.alphabetSize(); ++k) {
	  gsl_complex ck = gsl_complex_rect (0, 0);
	  for (AlphTok l = 0; l < model.alphabetSize(); ++l)
	    ck = gsl_complex_add
	      (ck,
	       gsl_complex_mul (eigenCounts[cpt][k][l],
				gsl_matrix_complex_get (evec[cpt], j, l)));
	  c = gsl_complex_add
	    (c,
	     gsl_complex_mul (gsl_matrix_complex_get (evecInv[cpt], k, i),
			      ck));
	}
	counts[i][j] = GSL_REAL(c) * (i == j ? 1 : gsl_matrix_get (model.subRate[cpt], i, j));
      }
    v.push_back (counts);
  }
  return v;
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

CachingRateModel::CachingRateModel (const RateModel& model, size_t precision, size_t flushSize)
  : RateModel (model),
    eigen (model),
    precision (precision),
    flushSize (flushSize)
{ }

string CachingRateModel::timeKey (double t) const {
  ostringstream o;
  o << scientific << setprecision(precision) << t;
  return o.str();
}

vguard<gsl_matrix*> CachingRateModel::getSubProbMatrix (double t) const {
  CachingRateModel* mutableThis = (CachingRateModel*) this;  // cast away const
  auto& mutableCache = mutableThis->cache;
  auto& mutableCount = mutableThis->count;
  const string k = timeKey (t);
  vguard<gsl_matrix*> m;
  if (cache.count (k)) {
    const auto& c = cache.at(k);
    for (int cpt = 0; cpt < components(); ++cpt)
      m.push_back (stl_to_gsl_matrix (c[cpt]));
  } else {
    m = eigen.getSubProbMatrix (t);
    if (mutableCount[k]++) {  // wait until the 2nd evaluation to start caching
      if (mutableCache.size() >= flushSize) {
	mutableCache.clear();
	mutableCount.clear();
      }
      vguard<vguard<vguard<double> > > c;
      for (int cpt = 0; cpt < components(); ++cpt)
	c.push_back (gsl_matrix_to_stl (m[cpt]));
      mutableCache[k] = c;
    }
  }
  return m;
}
