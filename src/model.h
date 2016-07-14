#ifndef MODEL_INCLUDED
#define MODEL_INCLUDED

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <string>

#include "gason.h"
#include "fastseq.h"
#include "logsumexp.h"
#include "alignpath.h"
#include "tree.h"

using namespace std;

#define DefaultDistanceMatrixIterations 100

#define DefaultCachingRateModelPrecision 5
#define DefaultCachingRateModelFlushSize 1000

struct AlphabetOwner {
  string alphabet;
  AlphabetOwner() { }
  AlphabetOwner (const string& a) { initAlphabet(a); }
  AlphabetOwner (const AlphabetOwner& ao) { initAlphabet(ao.alphabet); }
  void initAlphabet (const string& a);
  inline size_t alphabetSize() const { return alphabet.size(); }
  void readAlphabet (const JsonValue& json);
  UnvalidatedAlphTok tokenize (char c) const;
  AlphTok tokenizeOrDie (char c) const;
  bool isValidSymbol (char c) const;
  gsl_matrix* newAlphabetMatrix() const;
  gsl_vector* newAlphabetVector() const;

  void writeSubCounts (ostream& out, const vguard<vguard<double> >& rootCounts, const vguard<vguard<vguard<double> > >& subCountsAndWaitTimes, size_t indent = 0) const;
  void writeSubCountsComponent (ostream& out, const vguard<double>& rootCounts, const vguard<vguard<double> >& subCountsAndWaitTimes, size_t indent = 0) const;
};

struct RateModel : AlphabetOwner {
  double insRate, delRate, insExtProb, delExtProb;
  vguard<double> cptWeight;
  vguard<gsl_vector*> insProb;
  vguard<gsl_matrix*> subRate;

  RateModel();
  RateModel(const RateModel& model);
  ~RateModel();

  RateModel& operator= (const RateModel& model);

  int components() const { return cptWeight.size(); }

  void clear();
  void init (const string& alphabet);
  void read (const JsonValue& json);
  void write (ostream& out) const;
  void readComponent (const JsonMap& jm);
  void writeComponent (int cpt, ostream& out);

  vguard<gsl_vector*> getEqmProbVector() const;
  virtual vguard<gsl_matrix*> getSubProbMatrix (double t) const;

  double expectedSubstitutionRate() const;
  
  double mlDistance (const FastSeq& xGapped, const FastSeq& yGapped, int maxIterations = DefaultDistanceMatrixIterations) const;
  vguard<vguard<double> > distanceMatrix (const vguard<FastSeq>& gappedSeq, int maxIterations = DefaultDistanceMatrixIterations) const;
};

class EigenModel {
public:
  const RateModel& model;

  // indexed by component:
  vguard<gsl_vector_complex*> eval;
  vguard<gsl_matrix_complex*> evec;  // right eigenvectors
  vguard<gsl_matrix_complex*> evecInv;  // left eigenvectors

  EigenModel (const EigenModel& eigen);
  EigenModel (const RateModel& model);
  ~EigenModel();

  int components() const { return eval.size(); }

  vguard<gsl_matrix_complex*> eigenSubCount (double t) const;
  double getSubCount (int component, AlphTok a, AlphTok b, AlphTok i, AlphTok j, const gsl_matrix* sub, const gsl_matrix_complex* eSubCount) const;
  void accumSubCounts (int component, vguard<vguard<double> >& count, AlphTok a, AlphTok b, double weight, const gsl_matrix* sub, const gsl_matrix_complex* eSubCount) const;

  double getSubProb (int component, double t, AlphTok i, AlphTok j) const;
  vguard<gsl_matrix*> getSubProbMatrix (double t) const;
  gsl_matrix_complex* getRateMatrix (int component) const;
  gsl_matrix_complex* evecInv_evec (int component) const;

  vguard<vguard<vguard<double> > > getSubCounts (const vguard<vguard<vguard<gsl_complex> > >& eigenCounts) const;
  
private:
  vguard<vguard<gsl_complex> > ev, ev_t, exp_ev_t;

  vguard<bool> isReal;
  vguard<vguard<double> > realEval;
  vguard<vguard<vguard<double> > > realEvec, realEvecInv;
  vguard<vguard<double> > real_ev_t, real_exp_ev_t;

  void compute_exp_ev_t (double t, bool forceComplex = false);
  double getSubProbInner (int component, double t, AlphTok i, AlphTok j) const;
  
  EigenModel& operator= (const EigenModel&) = delete;
};

class CachingRateModel : public RateModel {
private:
  const size_t precision, flushSize;
  map<string,int> count;
  map<string,vector<vector<vector<double> > > > cache;
  string timeKey (double t) const;
  EigenModel eigen;
public:
  CachingRateModel (const RateModel& model, size_t precision = DefaultCachingRateModelPrecision, size_t flushSize = DefaultCachingRateModelFlushSize);
  vguard<gsl_matrix*> getSubProbMatrix (double t) const;
};

class ProbModel : public AlphabetOwner {
public:
  typedef enum { Start = 0,
		 Match = 0, Insert = 1, Delete = 2,
		 End = 3 } State;
  double t, ins, del, insExt, delExt;
  double insWait, delWait;
  vguard<double> cptWeight;
  vguard<gsl_vector*> insVec;
  vguard<gsl_matrix*> subMat;
  ProbModel (const RateModel& model, double t);
  ~ProbModel();
  int components() const { return cptWeight.size(); }
  double transProb (State src, State dest) const;
  void write (ostream& out) const;
  void writeComponent (int cpt, ostream& out) const;
  static State getState (bool parentUngapped, bool childUngapped);
private:
  ProbModel (const ProbModel&) = delete;
  ProbModel& operator= (const ProbModel&) = delete;
};

struct LogProbModel {
  typedef vguard<LogProb> LogProbVector;
  typedef vguard<vguard<LogProb> > LogProbMatrix;
  vguard<LogProb> logCptWeight;
  vguard<LogProbVector> logInsProb;
  vguard<LogProbMatrix> logSubProb;
  LogProbModel (const ProbModel& pm);
  int components() const { return logCptWeight.size(); }
};

struct IndelCounts {
  double ins, del, insExt, delExt, insTime, delTime;
  LogProb lp;
  IndelCounts (double pseudocount = 0, double pseudotime = 0);
  IndelCounts operator+ (const IndelCounts& c) const;
  IndelCounts operator* (double w) const;
  IndelCounts& operator+= (const IndelCounts& c);
  IndelCounts& operator*= (double w);
  void accumulateIndelCounts (const RateModel& model, double time, const AlignRowPath& parent, const AlignRowPath& child, double weight = 1.);
  void accumulateIndelCounts (const RateModel& model, const Tree& tree, const AlignPath& align, double weight = 1.);
  void writeJson (ostream& out, const size_t indent = 0) const;
  void read (const JsonValue& json);

  // helper to compute the expected wait time before an irreversible decay event,
  // conditioned on the event being known to have taken place
  static double decayWaitTime (double decayRate, double timeInterval);
};

struct EventCounts : AlphabetOwner {
  IndelCounts indelCounts;
  vguard<vguard<double> > rootCount;
  vguard<vguard<vguard<double> > > subCount;

  EventCounts() { }
  EventCounts (const AlphabetOwner& alph, int components, double pseudo = 0);

  int components() const { return rootCount.size(); }

  EventCounts operator+ (const EventCounts& c) const;
  EventCounts operator* (double w) const;
  EventCounts& operator+= (const EventCounts& c);
  EventCounts& operator*= (double w);

  void optimize (RateModel& model, bool fitIndelRates = true, bool fitSubstRates = true) const;
  
  void writeJson (ostream& out) const;
  void read (const JsonValue& json);
  void readComponent (const JsonMap& jm);

  double logPrior (const RateModel& model, bool includeIndelRates = true, bool includeSubstRates = true) const;
  double expectedLogLikelihood (const RateModel& model) const;
};

struct EigenCounts {
  IndelCounts indelCounts;
  vguard<vguard<double> > rootCount;
  vguard<vguard<vguard<gsl_complex> > > eigenCount;

  EigenCounts();
  EigenCounts (size_t components, size_t alphabetSize);

  int components() const { return rootCount.size(); }

  EigenCounts operator+ (const EigenCounts& c) const;
  EigenCounts operator* (double w) const;
  EigenCounts& operator+= (const EigenCounts& c);
  EigenCounts& operator*= (double w);

  void accumulateSubstitutionCounts (const RateModel& model, const Tree& tree, const vguard<FastSeq>& gapped, double weight = 1.);
  void accumulateCounts (const RateModel& model, const Alignment& align, const Tree& tree, bool updateIndelCounts = true, bool updateSubstCounts = true, double weight = 1.);

  EventCounts transform (const RateModel& model) const;
};

#endif /* MODEL_INCLUDED */
