#include "../src/kseq.h"
#include "../src/fastseq.h"
#include "../src/gason.h"
#include "../src/jsonutil.h"
#include "../src/knhx.h"
#include "../src/logger.h"
#include "../src/logsumexp.h"
#include "../src/vguard.h"
#include "../src/optparser.h"
#include "../src/recon.h"

// GNU --version
#define IDHIST_PROGNAME "idhist"
#define IDHIST_VERSION  "0.1"

struct ProgUsage : OptParser {
  ProgUsage (int argc, char** argv);
};

ProgUsage::ProgUsage (int argc, char** argv)
  : OptParser (argc, argv, IDHIST_PROGNAME, "{align,help,version} [options]")
{
  text = briefText
    + "\n"
    + "Commands:\n"
    + "\n"
    + " " + prog + " align seqs.fasta tree.newick  >alignment.fasta\n"
    + "\n";
}

int main (int argc, char** argv) {

  try {
    ProgUsage usage (argc, argv);
  
    const string command = usage.getCommand();

    if (command == "align") {

      Reconstructor recon;
      
      usage.implicitSwitches.push_back (string ("-seqs"));
      usage.implicitSwitches.push_back (string ("-tree"));
      
      deque<string>& argvec (usage.argvec);
      while (logger.parseLogArgs (argvec)
	     || recon.parseReconArgs (argvec)
	     || usage.parseUnknown())
	{ }

      Alignment align = recon.loadFilesAndReconstruct();
      writeFastaSeqs (cout, align.gapped());
      
    } else
      return usage.parseUnknownCommand (command, IDHIST_VERSION);

  } catch (...) {
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
