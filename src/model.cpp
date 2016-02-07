#include <math.h>

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_complex_math.h>

#include <algorithm>
#include <set>

#include "model.h"
#include "jsonutil.h"
#include "util.h"
#include "alignpath.h"
#include "logger.h"
#include "sumprod.h"

struct DistanceMatrixParams {
  const map<pair<AlphTok,AlphTok>,int>& pairCount;
  const RateModel& model;
  DistanceMatrixParams (const map<pair<AlphTok,AlphTok>,int>& pairCount, const RateModel& model)
    : pairCount(pairCount),
      model(model)
  { }
  double tML (int maxIterations) const;
  double negLogLike (double t) const;
};

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

RateModel::RateModel()
  : subRate(NULL),
    insProb(NULL)
{ }

RateModel::RateModel (const RateModel& model)
  : subRate(NULL),
    insProb(NULL)
{
  *this = model;
}

RateModel::~RateModel()
{
  if (subRate)
    gsl_matrix_free (subRate);
  if (insProb)
    gsl_vector_free (insProb);
}

RateModel& RateModel::operator= (const RateModel& model) {
  if (subRate)
    gsl_matrix_free (subRate);
  if (insProb)
    gsl_vector_free (insProb);

  ((AlphabetOwner&)*this) = model;
  insRate = model.insRate;
  delRate = model.delRate;
  insExtProb = model.insExtProb;
  delExtProb = model.delExtProb;
  insProb = model.newAlphabetVector();
  subRate = model.newAlphabetMatrix();
  CheckGsl (gsl_vector_memcpy (insProb, model.insProb));
  CheckGsl (gsl_matrix_memcpy (subRate, model.subRate));
  return *this;
}

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

  if (jm.contains("rootprob")) {
    const JsonMap insVec = jm.getObject ("rootprob");
    insProb = newAlphabetVector();
    gsl_vector_set_zero (insProb);

    for (auto ci : alphabet) {
      AlphTok i = (AlphTok) tokenize (ci);
      const string si (1, ci);
      if (insVec.contains (si))
	gsl_vector_set (insProb, i, insVec.getNumber(si));
    }
  } else
    insProb = getEqmProbVector();

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
  out << " \"rootprob\":" << endl;
  out << " {";
  for (AlphTok i = 0; i < alphabetSize(); ++i)
    out << (i>0 ? "," : "") << "\n  \"" << alphabet[i] << "\": " << gsl_vector_get(insProb,i);
  out << endl;
  out << " }," << endl;
  out << " \"subrate\":" << endl;
  out << " {" << endl;
  for (AlphTok i = 0; i < alphabetSize(); ++i) {
    out << "  \"" << alphabet[i] << "\": {";
    for (AlphTok j = 0; j < alphabetSize(); ++j)
      if (i != j)
	out << " \"" << alphabet[j] << "\": " << gsl_matrix_get(subRate,i,j) << (j < alphabetSize() - (i == alphabetSize() - 1 ? 2 : 1) ? "," : "");
    out << " }" << (i < alphabetSize() - 1 ? "," : "") << endl;
  }
  out << " }" << endl;
  out << "}" << endl;
}

gsl_vector* RateModel::getEqmProbVector() const {
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

gsl_matrix* RateModel::getSubProbMatrix (double t) const {
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
    insVec (model.newAlphabetVector()),
    subMat (model.getSubProbMatrix (t))
{
  CheckGsl (gsl_vector_memcpy (insVec, model.insProb));
}

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
  logInsProb = log_gsl_vector (pm.insVec);
  for (AlphTok i = 0; i < pm.alphabetSize(); ++i)
    for (AlphTok j = 0; j < pm.alphabetSize(); ++j)
      logSubProb[i][j] = log (gsl_matrix_get (pm.subMat, i, j));
}

double RateModel::mlDistance (const FastSeq& x, const FastSeq& y, int maxIterations) const {
  LogThisAt(4,"Estimating distance from " << x.name << " to " << y.name << endl);
  map<pair<AlphTok,AlphTok>,int> pairCount;
  Assert (x.length() == y.length(), "Sequences %s and %s have different lengths (%u, %u)", x.name.c_str(), y.name.c_str(), x.length(), y.length());
  for (size_t col = 0; col < x.length(); ++col) {
    const char ci = x.seq[col];
    const char cj = y.seq[col];
    if (!Alignment::isGap(ci) && !Alignment::isGap(cj))
      ++pairCount[pair<AlphTok,AlphTok> (tokenizeOrDie(ci), tokenizeOrDie(cj))];
  }
  if (LoggingThisAt(6)) {
    LogThisAt(6,"Counts:");
    for (const auto& pc : pairCount)
      LogThisAt(6," " << pc.second << "*" << alphabet[pc.first.first] << alphabet[pc.first.second]);
    LogThisAt(6,endl);
  }
  const DistanceMatrixParams dmp (pairCount, *this);
  return dmp.tML (maxIterations);
}

vguard<vguard<double> > RateModel::distanceMatrix (const vguard<FastSeq>& gappedSeq, int maxIterations) const {
  vguard<vguard<double> > dist (gappedSeq.size(), vguard<double> (gappedSeq.size()));
  for (size_t i = 0; i < gappedSeq.size() - 1; ++i)
    for (size_t j = i + 1; j < gappedSeq.size(); ++j)
      dist[i][j] = dist[j][i] = mlDistance (gappedSeq[i], gappedSeq[j]);
  if (LoggingThisAt(3)) {
    LogThisAt(3,"Distance matrix (" << dist.size() << " rows):" << endl);
    for (const auto& row : dist)
      LogThisAt(3,to_string_join(row) << endl);
  }
  return dist;
}

double distanceMatrixNegLogLike (double t, void *params) {
  const DistanceMatrixParams& dmp = *(const DistanceMatrixParams*) params;
  gsl_matrix* sub = dmp.model.getSubProbMatrix (t);
  double ll = 0;
  for (auto& pc : dmp.pairCount)
    ll += log (gsl_matrix_get(sub,pc.first.first,pc.first.second)) * (double) pc.second;
  gsl_matrix_free (sub);
  return -ll;
}

double DistanceMatrixParams::negLogLike (double t) const {
  return distanceMatrixNegLogLike (t, (void*) this);
}

double DistanceMatrixParams::tML (int maxIterations) const {
  const gsl_min_fminimizer_type *T;
  gsl_min_fminimizer *s;

  gsl_function F;
  F.function = &distanceMatrixNegLogLike;
  F.params = (void*) this;

  T = gsl_min_fminimizer_goldensection; // gsl_min_fminimizer_brent;
  s = gsl_min_fminimizer_alloc (T);

  double t;
  const double tLower = .0001, tUpper = 10 + tLower;
  const double llLower = negLogLike(tLower), llUpper = negLogLike(tUpper);
  LogThisAt(4,"tML: f(" << tLower << ") = " << llLower << ", f(" << tUpper << ") = " << llUpper << endl);

  const double maxSteps = 128;
  double nSteps;
  bool foundGuess = false;
  for (double step = 4; step > 1.0001 && !foundGuess; step = sqrt(step))
    for (double x = tLower + step; x < tUpper && !foundGuess; x += step) {
      const double ll = negLogLike(x);
      LogThisAt(4,"tML: f(" << x << ") = " << ll << endl);
      if (ll < llLower && ll < llUpper) {
	foundGuess = true;
	t = x;
      }
    }
  if (!foundGuess)
    return llLower < llUpper ? tLower : tUpper;

  LogThisAt(4,"Initializing with t=" << t << ", tLower=" << tLower << ", tUpper=" << tUpper << endl);
  gsl_min_fminimizer_set (s, &F, t, tLower, tUpper);

  for (int iter = 0; iter < maxIterations; ++iter) {
    (void) gsl_min_fminimizer_iterate (s);

    t = gsl_min_fminimizer_x_minimum (s);
    const double a = gsl_min_fminimizer_x_lower (s);
    const double b = gsl_min_fminimizer_x_upper (s);

    LogThisAt(5,"Iteration #" << iter+1 << " tML: f(" << t << ") = " << negLogLike(t) << endl);
    
    const int status = gsl_min_test_interval (a, b, 0.01, 0.0);  // converge to 1%
    if (status == GSL_SUCCESS)
      break;
  }

  gsl_min_fminimizer_free (s);

  return t;
}


void AlphabetOwner::writeSubCounts (ostream& out, const vguard<double>& rootCounts, const vguard<vguard<double> >& subCountsAndWaitTimes, size_t indent) const {
  const string ind (indent, ' ');
  out << ind << "{" << endl;
  out << ind << " \"root\":" << endl;
  out << ind << "  {";
  for (AlphTok i = 0; i < alphabetSize(); ++i)
    out << (i == 0 ? "" : ",") << endl << ind << "   \"" << alphabet[i] << "\": " << rootCounts[i];
  out << endl << ind << "  }," << endl;
  out << ind << " \"sub\":" << endl;
  out << ind << "  {";
  for (AlphTok i = 0; i < alphabetSize(); ++i) {
    out << (i == 0 ? "" : ",") << endl << ind << "   \"" << alphabet[i] << "\": {";
    for (AlphTok j = 0; j < alphabetSize(); ++j)
      if (i != j)
	out << (j == (i == 0 ? 1 : 0) ? "" : ",") << " " << "\"" << alphabet[j] << "\": " << subCountsAndWaitTimes[i][j];
    out << " }";
  }
  out << endl;
  out << ind << "  }," << endl;
  out << ind << " \"wait\":" << endl;
  out << ind << "  {";
  for (AlphTok i = 0; i < alphabetSize(); ++i)
    out << (i == 0 ? "" : ",") << endl << ind << "   \"" << alphabet[i] << "\": " << subCountsAndWaitTimes[i][i];
  out << endl;
  out << ind << "  }" << endl;
  out << ind << "}" << endl;
}

IndelCounts::IndelCounts()
  : ins(0), del(0), insExt(0), delExt(0), matchTime(0), delTime(0)
{ }

IndelCounts& IndelCounts::operator+= (const IndelCounts& c) {
  ins += c.ins;
  del += c.del;
  insExt += c.insExt;
  delExt += c.delExt;
  matchTime += c.matchTime;
  delTime += c.delTime;
  return *this;
}

IndelCounts& IndelCounts::operator*= (double w) {
  ins *= w;
  del *= w;
  insExt *= w;
  delExt *= w;
  matchTime *= w;
  delTime *= w;
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

EigenCounts::EigenCounts (size_t alphabetSize)
  : indelCounts(),
    rootCount (alphabetSize, 0),
    eigenCount (alphabetSize, vguard<gsl_complex> (alphabetSize, gsl_complex_rect(0,0)))
{ }

EigenCounts& EigenCounts::operator+= (const EigenCounts& c) {
  indelCounts += c.indelCounts;

  if (c.rootCount.size()) {
    if (rootCount.size())
      for (size_t n = 0; n < rootCount.size(); ++n)
	rootCount[n] += c.rootCount.at(n);
    else
      rootCount = c.rootCount;
  }

  if (c.eigenCount.size()) {
    if (eigenCount.size())
      for (size_t i = 0; i < eigenCount.size(); ++i)
	for (size_t j = 0; j < eigenCount[i].size(); ++j)
	  eigenCount[i][j] = gsl_complex_add (eigenCount[i][j], c.eigenCount.at(i).at(j));
    else
      eigenCount = c.eigenCount;
  }

  return *this;
}

EigenCounts& EigenCounts::operator*= (double w) {
  indelCounts *= w;
  for (auto& c : rootCount)
    c *= w;
  for (auto& sc : eigenCount)
    for (auto& c : sc)
      c = gsl_complex_mul_real (c, w);
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

EventCounts::EventCounts (const AlphabetOwner& alph)
  : AlphabetOwner (alph),
    indelCounts(),
    rootCount (alph.alphabetSize(), 0),
    subCount (alph.alphabetSize(), vguard<double> (alph.alphabetSize(), 0))
{ }

EventCounts& EventCounts::operator+= (const EventCounts& c) {
  Assert (alphabet == c.alphabet, "Alphabets don't match");

  indelCounts += c.indelCounts;

  for (size_t n = 0; n < rootCount.size(); ++n)
    rootCount[n] += c.rootCount.at(n);

  for (size_t i = 0; i < subCount.size(); ++i)
    for (size_t j = 0; j < subCount[i].size(); ++j)
      subCount[i][j] += c.subCount.at(i).at(j);

  return *this;
}

EventCounts& EventCounts::operator*= (double w) {
  indelCounts *= w;
  for (auto& c : rootCount)
    c *= w;
  for (auto& sc : subCount)
    for (auto& c : sc)
      c *= w;
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

void IndelCounts::accumulateIndelCounts (const AlignRowPath& parent, const AlignRowPath& child, double time, double weight) {
  enum { Match, Insert, Delete } state, next;
  state = Match;
  for (size_t col = 0; col < parent.size(); ++col) {
    if (parent[col] && child[col])
      next = Match;
    else if (parent[col])
      next = Delete;
    else if (child[col])
      next = Insert;
    else
      continue;
    switch (next) {
    case Match:
      matchTime += weight * time;
      break;
    case Insert:
      if (state == next)
	insExt += weight;
      else
	ins += weight;
      break;
    case Delete:
      if (state == next)
	delExt += weight;
      else {
	del += 1;
	delTime += weight * time;
      }
      break;
    default:
      break;
    }
    state = next;
  }
}

void IndelCounts::accumulateIndelCounts (const AlignPath& align, const Tree& tree, double weight) {
  for (TreeNodeIndex node = 0; node < tree.nodes() - 1; ++node)
    accumulateIndelCounts (align.at(tree.parentNode(node)), align.at(node), tree.branchLength(node), weight);
}

void EigenCounts::accumulateSubstitutionCounts (const RateModel& model, const Tree& tree, const vguard<FastSeq>& gapped, double weight) {
  AlignColSumProduct colSumProd (model, tree, gapped);

  EigenCounts c (model.alphabetSize());

  while (!colSumProd.alignmentDone()) {
    colSumProd.fillUp();
    colSumProd.fillDown();
    colSumProd.accumulateEigenCounts (c.rootCount, c.eigenCount);
    colSumProd.nextColumn();
  }

  c *= weight;
  *this += c;
}

void EigenCounts::accumulateCounts (const RateModel& model, const Alignment& align, const Tree& tree, double weight) {
  indelCounts.accumulateIndelCounts (align.path, tree, weight);
  accumulateSubstitutionCounts (model, tree, align.gapped(), weight);
}

EventCounts EigenCounts::transform (const RateModel& model) const {
  EventCounts c (model);
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
  out << "}" << endl;
}

void IndelCounts::writeJson (ostream& out, const size_t indent) const {
  const string ind (indent, ' ');
  out << ind << "{" << endl;
  out << ind << " \"ins\": " << ins << "," << endl;
  out << ind << " \"del\": " << del << "," << endl;
  out << ind << " \"insExt\": " << insExt << "," << endl;
  out << ind << " \"delExt\": " << delExt << "," << endl;
  out << ind << " \"matchTime\": " << matchTime << "," << endl;
  out << ind << " \"delTime\": " << delTime << endl;
  out << ind << "}";
}

void IndelCounts::read (const JsonValue& json) {
  JsonMap jm (json);
  ins = jm.getNumber("ins");
  del = jm.getNumber("del");
  insExt = jm.getNumber("insExt");
  delExt = jm.getNumber("delExt");
  matchTime = jm.getNumber("matchTime");
  delTime = jm.getNumber("delTime");
}

void EventCounts::read (const JsonValue& json) {
  AlphabetOwner alph;
  alph.readAlphabet (json);
  *this = EventCounts (alph);
  JsonMap jm (json);
  indelCounts.read (jm.getType("indel",JSON_OBJECT));
  JsonMap subBlock = jm.getObject("sub");
  JsonMap root = subBlock.getObject("root");
  JsonMap sub = subBlock.getObject("sub");
  JsonMap wait = subBlock.getObject("wait");
  for (AlphTok i = 0; i < alphabetSize(); ++i) {
    const string si (1, alphabet[i]);
    rootCount[i] = root.getNumber(si);
    subCount[i][i] = wait.getNumber(si);
    JsonMap sub_i = sub.getObject(si);
    for (AlphTok j = 0; j < alphabetSize(); ++j)
      if (i != j) {
	const string sj (1, alphabet[j]);
	subCount[i][j] = sub_i.getNumber(sj);
      }
  }
}
