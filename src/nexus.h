#ifndef NEXUS_INCLUDED
#define NEXUS_INCLUDED

#include "tree.h"
#include "fastseq.h"

// very crude Nexus format parser

#define DefaultNexusTreeName "tree"

struct NexusData {
  static const char gapChar, wildcardChar;

  vguard<FastSeq> gapped;
  vguard<string> rowName;
  string treeName;
  Tree tree;

  NexusData (const vguard<FastSeq>& matrix, const Tree& tree, bool convertGapsAndWildcards = true);
  NexusData (const string& nexusString);
  NexusData (istream& in);

  void read (const string& nexusString);
  void write (ostream& out) const;
};

#endif /* NEXUS_INCLUDED */
