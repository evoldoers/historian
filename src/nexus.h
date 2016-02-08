#ifndef NEXUS_INCLUDED
#define NEXUS_INCLUDED

#include "tree.h"
#include "fastseq.h"

// very crude Nexus format parser

struct NexusData {
  vguard<FastSeq> gapped;
  vguard<string> rowName;
  Tree tree;

  NexusData (const vguard<FastSeq>& matrix, const Tree& tree);
  NexusData (const string& nexusString);
  NexusData (istream& in);

  void read (const string& nexusString);
  void write (ostream& out) const;
};

#endif /* NEXUS_INCLUDED */
