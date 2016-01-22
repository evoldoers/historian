#ifndef MODEL_INCLUDED
#define MODEL_INCLUDED

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <string>

#include "gason.h"
#include "fastseq.h"

using namespace std;

struct AlphabetOwner {
  string alphabet;
  AlphabetOwner() { }
  AlphabetOwner (const string& a) : alphabet(a) { }
  AlphabetOwner (const AlphabetOwner& ao) : alphabet(ao.alphabet) { }
  inline size_t alphabetSize() const { return alphabet.size(); }
  void readAlphabet (const JsonValue& json);
  UnvalidatedAlphTok tokenize (char c) const;
  gsl_matrix* newAlphabetMatrix() const;
  gsl_vector* newAlphabetVector() const;
};

struct RateModel : AlphabetOwner {
  double insRate, delRate, insExtProb, delExtProb;
  gsl_matrix* subRate;

  RateModel();
  void read (const JsonValue& json);
  void write (ostream& out) const;
  gsl_vector* getEqmProb() const;
  gsl_matrix* getSubProb (double t) const;
};

struct ProbModel : AlphabetOwner {
  double t, ins, del, insExt, delExt;
  gsl_vector* insVec;
  gsl_matrix* subMat;
  ProbModel (const RateModel& model, double t);
  ~ProbModel();
  void write (ostream& out) const;
};

#endif /* MODEL_INCLUDED */
