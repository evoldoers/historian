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
#define HISTORIAN_PROGNAME "historian"
#define HISTORIAN_VERSION  "0.1"

struct ProgUsage : OptParser {
  ProgUsage (int argc, char** argv);
};

ProgUsage::ProgUsage (int argc, char** argv)
  : OptParser (argc, argv, HISTORIAN_PROGNAME, "{align,help,version} [options]")
{
  text = briefText
    + "\n"
    + "Commands:\n"
    + "\n"
    + " " + prog + " align seqs.fa [-tree tree.nh] [-model model.json] >alignment.fa\n"
    + "\n";
}

int main (int argc, char** argv) {

  ProgUsage usage (argc, argv);
  
  const string command = usage.getCommand();

  if (command == "align") {

    Reconstructor recon;

    usage.implicitSwitches.push_back (string ("-seqs"));

    deque<string>& argvec (usage.argvec);
    while (logger.parseLogArgs (argvec)
	   || recon.parseReconArgs (argvec)
	   || usage.parseUnknown())
      { }

    Alignment align = recon.loadFilesAndReconstruct();
    writeFastaSeqs (cout, align.gapped());
      
  } else
    return usage.parseUnknownCommand (command, HISTORIAN_VERSION);

  return EXIT_SUCCESS;
}
