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
    + "Reconstruction:\n"
    + "\n"
    + "  " + prog + " align seqs.fa [-tree tree.nh] [-model model.json] >alignment.fa\n"
    + "\n"
    + "File options:\n"
    + "   -seqs <file>    Specify unaligned sequences (FASTA)\n"
    + "   -guide <file>   Specify guide alignment (gapped FASTA)\n"
    + "   -tree <file>    Specify phylogeny (New Hampshire)\n"
    + "   -model <file>   Specify substitution & indel model (JSON)\n"
    + "\n"
    + "   -saveseqs <file>, -saveguide <file>, -savetree <file>, -savemodel <file>\n"
    + "                   Save various intermediate analysis results to files\n"
    + "\n"
    + "Reconstruction options:\n"
    + "   -band <n>       Size of band around guide alignment (default " + to_string(DefaultMaxDistanceFromGuide) + ")\n"
    + "   -noband         Turn off band, ignore guide alignment\n"
    + "   -minpost <p>    Posterior probability threshold for profile states (default " + TOSTRING(DefaultProfilePostProb) + ")\n"
    // Uncomment to show help for obsolescent random-sampling profile option:
    //    + "   -samples <n>    Sample profile states randomly\n"
    + "   -states <n>     Limit max number of states per profile\n"
    + "\n"
    // Uncomment to show help for obsolescent random-sampling profile option:
    //    + "Note -minpost and -samples select different profile construction strategies.\n"
    //    + " -minpost (the default) deterministically includes all states above a given\n"
    //    + "  posterior probability threshold, along with paths required to reach them.\n"
    //    + " -samples randomly samples state paths from the posterior distribution.\n"
    //    + "Both strategies are subject to the state limit imposed by -states.\n"
    //    + "\n"
    + "Guide alignment options:\n"
    + "   -kmatch <k>     Length of kmers for pre-filtering heuristic (default " + to_string(DEFAULT_KMER_LENGTH) + ")\n"
    + "   -kmatchn <n>    Threshold# of kmer matches to seed a diagonal\n"
    + "   -kmatchband <n> Size of DP band around kmer-matching diagonals (default " + to_string(DEFAULT_BAND_SIZE) + ")\n"
    + "   -kmatchmb <M>   Set kmer threshold to use M megabytes of memory\n"
    + "   -kmatchmax      Set kmer threshold to use all available memory (default)\n"
    + "   -kmatchoff      No kmer threshold, do full DP\n"
    + "\n"
    + "General options:\n"
    + "   -verbose, -vv, -vvv, -v4, -v5, etc.\n"
    // uncomment to document debug logging:
    //    + "   -log <function_name>\n"
    + "                   Various levels of logging (-nocolor for monochrome)\n"
    + "   -V, --version   Print GNU-style version info\n"
    + "   -h, --help      Print help message\n"
    + "   -seed <n>       Seed random number generator\n"
    + "\n"
    + "The method is that of phylogenetic transducers, as described in:\n"
    + " Westesson, Lunter, Paten & Holmes (2012).\n"
    + " Accurate Reconstruction of Insertion-Deletion Histories by\n"
    + " Statistical Phylogenetics.\n"
    + " PLoS One, DOI: 10.1371/journal.pone.0034572\n"
    + " http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0034572\n"
    + "\n"
    + "A longer, tutorial-style introduction is available here:\n"
    + " Phylogenetic automata, pruning, and multiple alignment\n"
    + " http://arxiv.org/abs/1103.4347"
    + "\n";
}

int main (int argc, char** argv) {

  ProgUsage usage (argc, argv);
  Reconstructor recon;

  const string command = usage.getCommand();
  deque<string>& argvec (usage.argvec);

  if (command == "align") {

    usage.implicitSwitches.push_back (string ("-seqs"));

    while (logger.parseLogArgs (argvec)
	   || recon.parseReconArgs (argvec)
	   || usage.parseUnknown())
      { }

    recon.loadReconFiles();
    recon.reconstruct();
    recon.writeRecon (cout);

    
  } else if (command == "count") {

    usage.implicitSwitches.push_back (string ("-recon"));
    usage.implicitSwitches.push_back (string ("-tree"));

    while (logger.parseLogArgs (argvec)
	   || recon.parseCountArgs (argvec)
	   || usage.parseUnknown())
      { }

    recon.loadCountFiles();
    recon.count();
    recon.writeCounts (cout);
    
  } else
    return usage.parseUnknownCommand (command, HISTORIAN_VERSION);

  
  return EXIT_SUCCESS;
}
