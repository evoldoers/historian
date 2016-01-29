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
    + "Reconstruction options:\n"
    + "   -seqs <file>    Specify unaligned sequences (FASTA)\n"
    + "   -guide <file>   Specify guide alignment (gapped FASTA)\n"
    + "   -tree <file>    Specify phylogeny (New Hampshire)\n"
    + "   -model <file>   Specify substitution & indel model (JSON)\n"
    + "\n"
    + "   -saveseqs <file>, -saveguide <file>, -savetree <file>, -savemodel <file>\n"
    + "                   Save intermediate analysis datasets to specified files\n"
    + "\n"
    + "   -band <n>       Size of band around guide alignment\n"
    + "   -noband         Turn off band, ignore guide alignment\n"
    + "   -samples <n>    Max number of alignments to sample when building profiles\n"
    + "   -states <n>     Max number of states allowed in any single profile\n"
    + "   -seed <n>       Seed random number generator\n"
    + "\n"
    + "Guide alignment options:\n"
    + "   -kmatch <k>     Length of kmers for pre-filtering heuristic (default " + to_string(DEFAULT_KMER_LENGTH) + ")\n"
    + "   -kmatchn <n>    Threshold# of kmer matches to seed a diagonal\n"
    + "   -kmatchband <n> Size of DP band around kmer-matching diagonals (default " + to_string(DEFAULT_BAND_SIZE) + ")\n"
    + "   -kmatchmb <M>   Set kmer threshold to use M megabytes of memory\n"
    + "   -kmatchmax      Set kmer threshold to use all available memory (default)\n"
    + "   -kmatchoff      No kmer threshold, do full DP\n"
    + "\n"
    + "The method is that of phylogenetic transducers, as described in\n"
    + "the following papers by Westesson, Lunter, Paten & Holmes:\n"
    + "\n"
    + "Accurate Reconstruction of Insertion-Deletion Histories by Statistical Phylogenetics\n"
    + "PLoS One, DOI: 10.1371/journal.pone.0034572\n"
    + "http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0034572\n"
    + "\n"
    + "Phylogenetic automata, pruning, and multiple alignment\n"
    + "http://arxiv.org/abs/1103.4347"
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
