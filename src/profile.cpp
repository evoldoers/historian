#include <sstream>
#include <math.h>
#include "profile.h"
#include "jsonutil.h"
#include "forward.h"
#include "alignpath.h"
#include "util.h"

#define WaitStateSuffix  ";"
#define ReadyStateSuffix "."

ProfileTransition::ProfileTransition()
  : lpTrans(-numeric_limits<double>::infinity())
{ }

ProfileState::ProfileState()
{ }

ProfileState::ProfileState (size_t components, AlphTok alphSize)
  : lpAbsorb (components, vguard<LogProb> (alphSize, -numeric_limits<double>::infinity()))
{ }

Profile::Profile (size_t components, const string& alphabet, const FastSeq& seq, AlignRowIndex rowIndex)
  : components (components),
    alphSize ((AlphTok) alphabet.size()),
    state (seq.length() + 2, ProfileState (components, (AlphTok) alphabet.size())),
    trans (seq.length() + 1),
    rootRowIndex (rowIndex)
{
  name = seq.name;
  state.front() = state.back() = ProfileState();  // start and end are null states
  state.front().name = "START";
  state.front().seqCoords[rowIndex] = 0;
  state.back().name = "END";
  state.back().seqCoords[rowIndex] = seq.length();
  set<char> invalidChars;
  int nInvalidToks = 0;
  for (size_t pos = 0; pos <= seq.seq.size(); ++pos) {
    ProfileTransition& t = trans[pos];
    t.src = pos;
    t.dest = pos + 1;
    t.lpTrans = 0;
    if (pos == seq.seq.size())
      state[pos].nullOut.push_back (pos);
    else
      state[pos].absorbOut.push_back (pos);
    state[pos+1].in.push_back (pos);
    if (pos < seq.seq.size()) {
      state[pos+1].name = string(1,seq.seq[pos]) + to_string(pos+1);
      state[pos+1].alignPath[rowIndex].push_back (true);
      state[pos+1].seqCoords[rowIndex] = pos + 1;
      auto& lpAbsorb = state[pos+1].lpAbsorb;
      for (auto& lpa: lpAbsorb)
	if (Alignment::isWildcard (seq.seq[pos]))
	  fill (lpa.begin(), lpa.end(), 0);
	else {
	  const UnvalidatedAlphTok tok = tokenize (seq.seq[pos], alphabet);
	  if (tok < 0) {
	    invalidChars.insert (seq.seq[pos]);
	    ++nInvalidToks;
	    fill (lpa.begin(), lpa.end(), 0);
	  } else
	    lpa[tok] = 0;
	}
    }
  }
  this->seq[rowIndex] = seq.seq;

  if (nInvalidToks)
    Warn("%s (%s) found in sequence %s", plural(nInvalidToks,"invalid character").c_str(), join(invalidChars).c_str(), seq.name.c_str());
  
  assertTransitionsConsistent();
  assertSeqCoordsConsistent();
  assertAllStatesWaitOrReady();
  assertPathToEndExists();
}

Profile Profile::leftMultiply (const vguard<gsl_matrix*>& sub) const {
  Profile prof (*this);
  for (ProfileStateIndex i = 0; i < size(); ++i)
    if (!state[i].isNull())
      for (int cpt = 0; cpt < components; ++cpt) {
	for (AlphTok c = 0; c < alphSize; ++c) {
	  LogProb lp = -numeric_limits<double>::infinity();
	  for (AlphTok d = 0; d < alphSize; ++d)
	    lp = log_sum_exp (lp, log (gsl_matrix_get(sub[cpt],c,d)) + state[i].lpAbsorb[cpt][d]);
	  prof.state[i].lpAbsorb[cpt][c] = lp;
	}
      }
  return prof;
}

const ProfileTransition* Profile::getTrans (ProfileStateIndex src, ProfileStateIndex dest) const {
  for (auto t : state[dest].in)
    if (trans[t].src == src)
      return &trans[t];
  return NULL;
}

map<AlignRowIndex,char> Profile::alignColumn (ProfileStateIndex s) const {
  map<AlignRowIndex,char> col;
  for (auto& row_path : state[s].alignPath)
    if (row_path.second.size() && row_path.second.front()) {
      if (state[s].seqCoords.count(row_path.first))
	col[row_path.first] = seq.at(row_path.first).at(state[s].seqCoords.at(row_path.first) - 1);
      else
	col[row_path.first] = Alignment::wildcardChar;
    }
  return col;
}

LogProb Profile::calcSumPathAbsorbProbs (const vguard<LogProb>& logCptWeight, const vguard<vguard<LogProb> >& logInsProb, const char* tag) {
  vguard<LogProb> lpCumAbs (state.size(), -numeric_limits<double>::infinity());
  lpCumAbs[0] = 0;
  for (ProfileStateIndex pos = 1; pos < state.size(); ++pos) {
    LogProb lpAbs = 0;
    if (!state[pos].isNull()) {
      lpAbs = -numeric_limits<double>::infinity();
      for (int cpt = 0; cpt < components; ++cpt)
	log_accum_exp (lpAbs, logCptWeight[cpt] + logInnerProduct (logInsProb[cpt], state[pos].lpAbsorb[cpt]));
    }
    for (auto ti : state[pos].in) {
      const ProfileTransition& t = trans[ti];
      Assert (t.src < pos, "Transition #%u from %u -> %u is not toposorted", ti, t.src, t.dest);
      log_accum_exp (lpCumAbs[pos], lpCumAbs[t.src] + t.lpTrans + lpAbs);
    }
    if (tag != NULL)
      state[pos].meta[string(tag)] = to_string(lpCumAbs[pos]);
  }
  return lpCumAbs.back();
}

string alignPathJson (const AlignPath& a) {
  string s = "[";
  for (auto& row_path : a) {
    if (s.size() > 1)
      s += ",";
    s += " [ " + to_string(row_path.first) + ", \"";
    for (auto col : row_path.second)
      s += (col ? Alignment::wildcardChar : Alignment::gapChar);
    s += "\" ]";
  }
  s += " ]";
  return s;
}

void Profile::writeJson (ostream& out) const {
  out << "{" << endl;
  if (name.size())
    out << " \"name\": \"" << name << "\"," << endl;
  if (meta.size())
    out << " \"meta\": " << JsonUtil::toString(meta,2) << "," << endl;
  out << " \"alphSize\": " << alphSize << "," << endl;
  out << " \"state\": [" << endl;
  for (ProfileStateIndex s = 0; s < state.size(); ++s) {
    out << "  {" << endl;
    out << "   \"n\": " << s << "," << endl;
    if (state[s].name.size())
      out << "   \"name\": \"" << state[s].name << "\"," << endl;
    if (state[s].meta.size())
      out << "   \"meta\": " << JsonUtil::toString(state[s].meta,4) << "," << endl;
    if (state[s].alignPath.size())
      out << "   \"path\": " << alignPathJson(state[s].alignPath) << "," << endl;
    if (state[s].seqCoords.size()) {
      out << "   \"seqPos\": [";
      size_t nSeqPos = 0;
      for (const auto& s_c : state[s].seqCoords)
	out << (nSeqPos++ ? ", " : " ") << "[ " << s_c.first << ", " << s_c.second << " ]";
      out << " ]," << endl;
    }
    if (!state[s].isNull()) {
      out << "   \"lpAbsorb\": [";
      for (int cpt = 0; cpt < components; ++cpt) {
	out << (cpt > 0 ? ", " : "") << "[";
	for (AlphTok a = 0; a < alphSize; ++a)
	  out << (a > 0 ? ", " : " ") << JsonUtil::toString (state[s].lpAbsorb[cpt][a]);
	out << " ]";
      }
      out << "]," << endl;
    }
    out << "   \"trans\": [";
    set<ProfileTransitionIndex> s_out (state[s].nullOut.begin(), state[s].nullOut.end());
    s_out.insert (state[s].absorbOut.begin(), state[s].absorbOut.end());
    bool first_t = true;
    for (auto ti : s_out) {
      const ProfileTransition& tr = trans[ti];
      if (!first_t)
	out << ",\n             ";     
      first_t = false;
      out << " { \"to\": " << tr.dest << ",";
      out << " \"lpTrans\": " << JsonUtil::toString (tr.lpTrans);
      if (tr.alignPath.size())
	out << ", \"path\": " << alignPathJson(tr.alignPath);
      out << " }";
    }
    out << " ]" << endl;
    out << "  }";
    if (s < state.size() - 1)
      out << ",";
    out << endl;
  }
  out << " ]" << endl;
  out << "}" << endl;
}

string Profile::toJson() const {
  ostringstream s;
  writeJson (s);
  return s.str();
}

string trans2state (const vguard<ProfileTransition>& trans, const vguard<ProfileTransitionIndex>& idx, bool useSrc) {
  vguard<ProfileStateIndex> v;
  for (auto i: idx)
    v.push_back (useSrc ? trans[i].src : trans[i].dest);
  return string("(") + to_string_join(v,",") + ")";
}

string Profile::tinyDescription (ProfileStateIndex i) const {
  ostringstream out;
  const ProfileState& s = state[i];
  out << "("
      << (s.isStart() ? "Start" : (s.isNull() ? "Null" : "Emit"))
      << ","
      << (s.isReady() ? "Ready" : (s.isWait() ? "Wait" : "Irregular"));
  if (!s.in.empty())
    out << ",in:" << trans2state (trans, s.in, true);
  if (!s.nullOut.empty())
    out << ",nullOut:" << trans2state (trans, s.nullOut, false);
  if (!s.absorbOut.empty())
    out << ",absorbOut:" << trans2state (trans, s.absorbOut, false);
  out << ")";
  return out.str();
}

void Profile::assertSeqCoordsConsistent() const {
  for (const auto& t: trans)
    ProfileState::assertSeqCoordsConsistent (state[t.src].seqCoords, state[t.dest], t.alignPath);
}

void ProfileState::assertSeqCoordsConsistent (const SeqCoords& srcCoords, const ProfileState& dest, const AlignPath& transPath) {
  assertSeqCoordsConsistent (srcCoords, dest.seqCoords, transPath, dest.alignPath);
}

void ProfileState::assertSeqCoordsConsistent (const SeqCoords& srcCoords, const SeqCoords& destCoords, const AlignPath& transPath, const AlignPath& destPath) {
  SeqCoords seqCoords = srcCoords;
    for (const auto& rp: transPath)
      seqCoords[rp.first] += alignPathResiduesInRow (rp.second);
    for (const auto& rp: destPath)
      seqCoords[rp.first] += alignPathResiduesInRow (rp.second);
    for (const auto& sc: destCoords) {
      Assert (seqCoords.count(sc.first), "Missing coordinate for sequence %d", sc.first);
      Assert (seqCoords.at(sc.first) == sc.second,
	      "Sequence coord %d: source state (%d) + transition path (%d) + dest state path (%d) != dest state (%d)",
	      sc.first,
	      srcCoords.count(sc.first) ? srcCoords.at(sc.first) : 0,
	      transPath.count(sc.first) ? alignPathResiduesInRow(transPath.at(sc.first)) : 0,
	      destPath.count(sc.first) ? alignPathResiduesInRow(destPath.at(sc.first)) : 0,
	      sc.second);
    }
}

void Profile::assertAllStatesWaitOrReady() const {
  for (auto& s: state)
    Assert (s.isReady() || s.isWait(), "State %s has %d null transitions and %d absorbing transitions, so is neither Wait nor Ready", s.name.c_str(), s.nullOut.size(), s.absorbOut.size());
}

Profile Profile::addReadyStates() const {
  vguard<ProfileStateIndex> old2newStateIndex (size());
  Profile prof;
  prof.alphSize = alphSize;
  prof.components = components;
  prof.name = name;
  prof.meta = meta;
  prof.seq = seq;
  prof.trans = trans;
  prof.rootRowIndex = rootRowIndex;
  vguard<ProfileState> profState (state);
  for (ProfileStateIndex s = 0, n = 0; s < size(); ++s) {
    old2newStateIndex[s] = n++;
    if (!state[s].isReady() && !state[s].isWait()) {
      ProfileState readyState;
      ProfileTransition readyTrans;
      const ProfileStateIndex oldReadyStateIdx = profState.size();
      const ProfileStateIndex newReadyStateIdx = n++;
      const ProfileTransitionIndex readyTransIdx = prof.trans.size();
      profState[s].name += WaitStateSuffix;
      readyState.name = state[s].name + ReadyStateSuffix;
      readyState.meta = state[s].meta;
      readyState.seqCoords = state[s].seqCoords;
      swap (profState[s].absorbOut, readyState.absorbOut);
      for (auto t: readyState.absorbOut)
	prof.trans[t].src = oldReadyStateIdx;
      readyTrans.src = s;
      readyTrans.dest = oldReadyStateIdx;
      readyTrans.lpTrans = 0;
      profState[s].nullOut.push_back (readyTransIdx);
      readyState.in.push_back (readyTransIdx);
      profState.push_back (readyState);
      prof.trans.push_back (readyTrans);
      old2newStateIndex.push_back (newReadyStateIdx);
    }
  }
  prof.state = vguard<ProfileState> (profState.size());
  for (ProfileStateIndex s = 0; s < (int) profState.size(); ++s)
    swap (profState[s], prof.state[old2newStateIndex[s]]);
  for (auto& t: prof.trans) {
    t.src = old2newStateIndex[t.src];
    t.dest = old2newStateIndex[t.dest];
  }
  for (const auto& ss: equivAbsorbState)
    prof.equivAbsorbState[old2newStateIndex[ss.first]] = old2newStateIndex[ss.second];

  prof.assertTransitionsConsistent();
  prof.assertAllStatesWaitOrReady();
  prof.assertPathToEndExists();

  return prof;
}

void Profile::assertTransitionsConsistent() const {
  for (ProfileStateIndex i = 0; i < (int) state.size(); ++i) {
    const ProfileState& s = state[i];
    for (ProfileTransitionIndex t: s.in)
      Assert (trans[t].dest == i, "Incoming transition destination index doesn't match state index");
    for (ProfileTransitionIndex t: s.nullOut)
      Assert (trans[t].src == i, "Null transition source index doesn't match state index");
    for (ProfileTransitionIndex t: s.absorbOut)
      Assert (trans[t].src == i, "Absorbing transition source index doesn't match state index");
  }
}

void Profile::assertPathToEndExists() const {
  (void) examplePathToEnd();
}

vguard<ProfileStateIndex> Profile::examplePathToEnd() const {
  vguard<bool> fromStart (state.size(), false);
  vguard<ProfileStateIndex> prev (state.size(), 0);
  fromStart.front() = true;
  for (ProfileStateIndex i = 0; i < (int) state.size(); ++i)
    if (fromStart[i]) {
      const ProfileState& s = state[i];
      for (ProfileTransitionIndex t: s.nullOut) {
	Assert (trans[t].dest > i, "Null transition violates toposort");
	fromStart[trans[t].dest] = true;
	prev[trans[t].dest] = i;
      }
      for (ProfileTransitionIndex t: s.absorbOut) {
	Assert (trans[t].dest > i, "Absorbing transition violates toposort");
	fromStart[trans[t].dest] = true;
	prev[trans[t].dest] = i;
      }
    }
  Assert (fromStart.back(), "No path from start to end");
  vguard<ProfileStateIndex> revPath;
  for (ProfileStateIndex j = state.size() - 1; j != 0; j = prev[j])
    revPath.push_back (j);
  revPath.push_back (0);
  return vguard<ProfileStateIndex> (revPath.rbegin(), revPath.rend());
}

bool Profile::isEmpty() const {
  for (auto s: state)
    if (!s.isNull())
      return false;
  return true;
}
