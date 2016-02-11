#include <algorithm>
#include <iomanip>
#include "nexus.h"
#include "jsonutil.h"
#include "util.h"
#include "alignpath.h"

const char NexusData::gapChar = '-';
const char NexusData::wildcardChar = '?';

NexusData::NexusData (const vguard<FastSeq>& matrix, const Tree& tree)
  : gapped (matrix),
    tree (tree),
    treeName (DefaultNexusTreeName)
{ }

void NexusData::convertNexusToAlignment() {
  for (auto& fs : gapped)
    for (auto& c : fs.seq)
      if (c == gapChar)
	c = Alignment::gapChar;
      else if (c == wildcardChar)
	c = Alignment::wildcardChar;
}

void NexusData::convertAlignmentToNexus() {
  for (auto& fs : gapped)
    for (auto& c : fs.seq)
      if (c == Alignment::gapChar)
	c = gapChar;
      else if (c == Alignment::wildcardChar)
	c = wildcardChar;
}

NexusData::NexusData (const string& nexusString) {
  read (nexusString);
}

NexusData::NexusData (istream& in) {
  read (JsonUtil::readStringFromStream (in, true));
}

void NexusData::read (const string& nexusString) {
  string preproc;
  preproc.reserve (nexusString.size());
  enum { NoComment, HashComment, SquareBracketComment } state = NoComment;
  for (auto c : nexusString) {
    switch (state) {
    case NoComment:
      if (c == '#')
	state = HashComment;
      else if (c == '[')
	state = SquareBracketComment;
      else
	preproc.push_back (c);
      break;
    case HashComment:
      if (c == '\n')
	state = NoComment;
      break;
    case SquareBracketComment:
      if (c == ']')
	state = NoComment;
      break;
    default:
      break;
    }
  }
  
  const vguard<string> line = split (preproc, ";");
  enum { NoBlock, Data, Tree } block = NoBlock;
  map<string,string> seq;
  vguard<string> rowName;
  for (size_t l = 0; l < line.size(); ++l) {
    vguard<string> tok = split (line[l]);
    if (tok.size()) {
      const string cmd = toupper(tok[0]);
      switch (block) {
      case NoBlock:
	if (tok.size() == 2 && cmd == "BEGIN") {
	  if (toupper(tok[1]) == "DATA") {
	    Require (seq.empty(), "Multiple DATA blocks in Nexus file");
	    block = Data;
	  } else if (toupper(tok[1]) == "TREE" || toupper(tok[1]) == "TREES")
	    block = Tree;
	}
	break;

      case Data:
	if (tok.size() == 1 && cmd == "END")
	  block = NoBlock;
	else if (cmd == "MATRIX") {
	  Require (tok.size() % 2 == 1, "MATRIX block does not have an even number of fields");
	  for (size_t n = 1; n < tok.size(); n += 2) {
	    if (seq.find (tok[n]) == seq.end())
	      rowName.push_back (tok[n]);
	    seq[tok[n]] += tok[n+1];
	  }
	}
	break;

      case Tree:
	if (tok.size() == 1 && cmd == "END")
	  block = NoBlock;
	else if (cmd == "TREE" && tok.size() == 4 && tok[2] == "=") {
	  Require (treeName.empty(), "Multiple trees in Nexus file");
	  treeName = tok[1];
	  tree.parse (tok[3] + ";");
	}
	break;

      default:
	break;
      }
    }
  }
  Require (rowName.size() > 0, "No sequence data found in Nexus file");
  Require (tree.nodes() > 0, "No tree found in Nexus file");
  for (size_t n = 0; n < rowName.size(); ++n) {
    FastSeq fs;
    fs.name = rowName[n];
    fs.seq = seq[rowName[n]];
    gapped.push_back (fs);
  }
}

void NexusData::write (ostream& out) const {
  out << "#NEXUS" << endl;
  out << "BEGIN DATA;" << endl;
  if (gapped.size()) {
    out << "DIMENSIONS NTAX=" << gapped.size() << " NCHAR=" << gapped[0].length() << ";" << endl;
    out << "MATRIX" << endl;
    size_t w = 0;
    for (const auto& fs : gapped)
      w = max (w, fs.name.size());
    for (const auto& fs : gapped)
      out << setw(w+1) << left << fs.name << fs.seq << endl;
    out << ";" << endl;
  }
  out << "END;" << endl;
  out << "BEGIN TREES;" << endl;
  out << "TREE " << treeName << " = " << tree.toString() << endl;
  out << "END;" << endl;
}
