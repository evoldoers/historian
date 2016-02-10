#include <algorithm>
#include <iomanip>
#include "stockholm.h"
#include "regexmacros.h"
#include "util.h"

// POSIX basic regular expressions
const regex seq_re (RE_WHITE_OR_EMPTY RE_GROUP(RE_PLUS(RE_NONWHITE_CHAR_CLASS)) RE_WHITE_NONEMPTY RE_GROUP(RE_PLUS(RE_NONWHITE_CHAR_CLASS)) RE_WHITE_OR_EMPTY, regex_constants::basic);
const regex gf_re (RE_WHITE_OR_EMPTY "#=GF" RE_WHITE_NONEMPTY RE_GROUP(RE_PLUS(RE_NONWHITE_CHAR_CLASS)) RE_WHITE_NONEMPTY RE_GROUP(RE_DOT_STAR), regex_constants::basic);
const regex gc_re (RE_WHITE_OR_EMPTY "#=GC" RE_WHITE_NONEMPTY RE_GROUP(RE_PLUS(RE_NONWHITE_CHAR_CLASS)) RE_WHITE_NONEMPTY RE_GROUP(RE_DOT_STAR), regex_constants::basic);
const regex gr_re (RE_WHITE_OR_EMPTY "#=GR" RE_WHITE_NONEMPTY RE_GROUP(RE_PLUS(RE_NONWHITE_CHAR_CLASS)) RE_WHITE_NONEMPTY RE_GROUP(RE_PLUS(RE_NONWHITE_CHAR_CLASS)) RE_WHITE_NONEMPTY RE_GROUP(RE_DOT_STAR), regex_constants::basic);
const regex gs_re (RE_WHITE_OR_EMPTY "#=GS" RE_WHITE_NONEMPTY RE_GROUP(RE_PLUS(RE_NONWHITE_CHAR_CLASS)) RE_WHITE_NONEMPTY RE_GROUP(RE_PLUS(RE_NONWHITE_CHAR_CLASS)) RE_WHITE_NONEMPTY RE_GROUP(RE_DOT_STAR), regex_constants::basic);
const regex hash_re (RE_WHITE_OR_EMPTY "#" RE_DOT_STAR, regex_constants::basic);
const regex divider_re (RE_WHITE_OR_EMPTY "//" RE_WHITE_OR_EMPTY, regex_constants::basic);

void Stockholm::read (istream& in) {
  gf.clear();
  gc.clear();
  gs.clear();
  gr.clear();
  gapped.clear();

  smatch sm;
  map<string,string> seq;
  vguard<string> rowName;
  while (in && !in.eof()) {
    string line;
    getline(in,line);
    if (regex_match (line, sm, seq_re)) {
      if (!seq.count (sm.str(1)))
	rowName.push_back (sm.str(1));
      seq[sm.str(1)] += sm.str(2);
    } else if (regex_match (line, sm, gf_re))
      gf[sm.str(1)].push_back (sm.str(2));
    else if (regex_match (line, sm, gc_re))
      gc[sm.str(1)] += sm.str(2);
    else if (regex_match (line, sm, gr_re))
      gr[sm.str(2)][sm.str(1)] += sm.str(3);
    else if (regex_match (line, sm, gs_re))
      gs[sm.str(2)][sm.str(1)].push_back (sm.str(3));
    else if (regex_match (line, hash_re))
      continue;
    else if (regex_match (line, divider_re))
      break;
    else
      Warn ("Unrecognized line in Stockholm file: %s", line.c_str());
  }
  for (const auto& name : rowName) {
    FastSeq fs;
    fs.name = name;
    fs.seq = seq[name];
    gapped.push_back (fs);
  }
}

void Stockholm::write (ostream& out, size_t charsPerRow) const {
  int nw = 0, tw = 0, cols = columns();
  set<string> names;
  for (auto& fs : gapped) {
    nw = max (nw, (int) fs.name.size());
    names.insert (fs.name);
  }
  for (auto& tag_gf : gf)
    tw = max (tw, (int) tag_gf.first.size());
  for (auto& tag_gc : gc) {
    tw = max (tw, (int) tag_gc.first.size());
    cols = max (cols, (int) tag_gc.second.size());
  }
  for (auto& tag_gs : gs) {
    tw = max (tw, (int) tag_gs.first.size());
    for (auto& name_gs : tag_gs.second)
      nw = max (nw, (int) name_gs.first.size());
  }
  for (auto& tag_gr : gr) {
    tw = max (tw, (int) tag_gr.first.size());
    for (auto& name_gr : tag_gr.second) {
      nw = max (nw, (int) name_gr.first.size());
      cols = max (cols, (int) name_gr.second.size());
    }
  }
  const int w = tw > 0 ? (nw + tw + 6) : nw;
  
  out << "# STOCKHOLM 1.0" << endl;
  for (auto& tag_gf : gf)
    for (auto& line : tag_gf.second)
      out << setw(tw+6) << left << (string("#=GF ") + tag_gf.first) << line << endl;
  for (auto& tag_gs : gs) {
    for (auto& fs : gapped)
      if (tag_gs.second.count (fs.name))
	for (auto& line : tag_gs.second.at(fs.name))
	  out << "#=GS" << right << setw(nw+1) << fs.name << " " << left << setw(tw+1) << tag_gs.first << line << endl;
    for (auto& name_gs : tag_gs.second)
      if (!names.count (name_gs.first))
	for (auto& line : name_gs.second)
	  out << "#=GS" << right << setw(nw+1) << name_gs.first << " " << left << setw(tw+1) << tag_gs.first << line << endl;
  }
  const int colStep = min (MinStockholmCharsPerRow, DefaultStockholmRowLength - w - 1);
  for (int col = 0; col < cols; col += colStep) {
    for (auto& tag_gc : gc)
      out << "#=GC" << right << setw(w-4) << tag_gc.first << " " << tag_gc.second.substr(col,colStep) << endl;
    for (auto& fs : gapped) {
      out << left << setw(w+1) << fs.name << fs.seq.substr(col,colStep) << endl;
      for (auto& tag_gr : gr)
	if (tag_gr.second.count (fs.name))
	  out << "#=GR" << right << setw(nw+1) << fs.name << " " << left << setw(tw+1) << tag_gr.first << tag_gr.second.at(fs.name).substr(col,colStep) << endl;
    }
    for (auto& tag_gr : gr)
      for (auto& name_gr : tag_gr.second)
	if (!names.count (name_gr.first))
	  out << "#=GR" << right << setw(nw+1) << name_gr.first << " " << left << setw(tw+1) << tag_gr.first << tag_gr.second.at(name_gr.first).substr(col,colStep) << endl;
  }
  out << "//" << endl;
}

void Stockholm::setTree (const Tree& tree, const char* tag) {
  gf[string(tag)].push_back (tree.toString());
}

size_t Stockholm::rows() const {
  return gapped.size();
}

size_t Stockholm::columns() const {
  return gappedSeqColumns (gapped);
}

AlignPath Stockholm::path() const {
  Alignment a (gapped);
  return a.path;
}
