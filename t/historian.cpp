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
  : OptParser (argc, argv, HISTORIAN_PROGNAME, "{recon[struct],count,fit,help,version} [options]")
{
  text = briefText
    + "\n"
    + "EXAMPLES\n"
    + "\n"
    + "Reconstruction:\n"
    + "  " + prog + " recon seqs.fa [-tree tree.nh] -output fasta >reconstruction.fa\n"
    + "  " + prog + " recon -guide guide.fa [-tree tree.nh] >reconstruction.stk\n"
    + "  " + prog + " recon guide.stk >reconstruction.stk\n"
    + "  " + prog + " recon data.nex -output nexus >reconstruction.nex\n"
    + "\n"
    + "Event counting:\n"
    + "  " + prog + " count seqs.fa [-tree tree.nh] [-model model.json] >counts.json\n"
    + "  " + prog + " count -guide guide.fa [-tree tree.nh] >counts.json\n"
    + "  " + prog + " count -recon reconstruction.fa -tree tree.nh >counts.json\n"
    + "\n"
    + "Model fitting:\n"
    + "  " + prog + " fit seqs.fa >newmodel.json\n"
    + "  " + prog + " fit -counts counts.json >newmodel.json\n"
    + "\n"
    + "Commands can be abbreviated to single letters, like so:\n"
    + "  " + prog + " r seqs.fa >reconstruction.stk\n"
    + "  " + prog + " c seqs.fa >counts.json\n"
    + "  " + prog + " f -counts counts.json >model.json\n"
    + "(etc.)\n"
    + "\n"
    + "OPTIONS\n"
    + "\n"
    + "Reconstruction file I/O options:\n"
    + "  -auto <file>    Auto-detect file format and purpose\n"
    + "  -model <file>   Specify substitution & indel model (JSON)\n"
    + "  -seqs <file>    Specify unaligned sequences (FASTA)\n"
    + "  -guide <file>   Specify guide alignment (gapped FASTA)\n"
    + "  -tree <file>    Specify phylogeny (New Hampshire)\n"
    + "  -nexus <file>, -stockholm <file>\n"
    + "                  Specify phylogeny & guide alignment together\n"
    + "\n"
    + "  -saveguide <f>  Save guide alignment to file\n"
    + "  -output (nexus|fasta|stockholm)\n"
    + "                  Specify output format (default is Stockholm)\n"
    + "\n"
    + "Reconstruction algorithm options:\n"
    + "  -ancseq         Predict ancestral sequences (default is to leave them as " + Alignment::wildcardChar + "'s)\n"
    + "  -band <n>       Size of band around guide alignment (default " + to_string(DefaultMaxDistanceFromGuide) + ")\n"
    + "  -noband         Turn off band, ignore guide alignment\n"
    + "  -minpost <p>    Posterior prob. threshold for profile states (default " + TOSTRING(DefaultProfilePostProb) + ")\n"
    // Uncomment to show help for obsolescent random-sampling profile option:
    //    + "  -samples <n>    Sample profile states randomly\n"
    + "  -states <n>     Limit max number of states per profile\n"
    + "\n"
    // Uncomment to show help for obsolescent random-sampling profile option:
    //    + "Note -minpost and -samples select different profile construction strategies.\n"
    //    + " -minpost (the default) deterministically includes all states above a given\n"
    //    + "  posterior probability threshold, along with paths required to reach them.\n"
    //    + " -samples randomly samples state paths from the posterior distribution.\n"
    //    + "Both strategies are subject to the state limit imposed by -states.\n"
    //    + "\n"
    + "Guide alignment & tree estimation options:\n"
    + "  -allvsall       Try all pairwise alignments, not just a random spanning graph\n"
    + "  -kmatch <k>     Length of kmers for pre-filtering heuristic (default " + to_string(DEFAULT_KMER_LENGTH) + ")\n"
    + "  -kmatchn <n>    Threshold# of kmer matches to seed a diagonal\n"
    + "  -kmatchband <n> Size of DP band around kmer-matching diagonals (default " + to_string(DEFAULT_BAND_SIZE) + ")\n"
    + "  -kmatchmb <M>   Set kmer threshold to use M megabytes of memory\n"
    + "  -kmatchmax      Set kmer threshold to use all available memory (default)\n"
    + "  -kmatchoff      No kmer threshold, do full DP\n"
    + "  -upgma          Use UPGMA, not neighbor-joining, to estimate tree\n"
    + "\n"
    + "Model-fitting and event-counting options:\n"
    + "  -recon <file>, -nexusrecon <file>, -stockrecon <file>\n"
    + "                  Use precomputed reconstruction (FASTA/NEXUS/Stockholm)\n"
    + "  -mininc <n>     EM convergence threshold as relative log-likelihood increase\n"
    + "                   (default is " + TOSTRING(DefaultMinEMImprovement) + ")\n"
    + "  -maxiter <n>    Max number of EM iterations (default " + to_string(DefaultMaxEMIterations) + ")\n"
    + "  -nolaplace      Do not add Laplace +1 pseudocounts during model-fitting\n"
    + "  -fixsubrates    Do not estimate substitution rates or initial composition\n"
    + "  -fixgaprates    Do not estimate indel rates or length distributions\n"
    + "\n"
    + "General options:\n"
    + "  -verbose, -vv, -vvv, -v4, -v5, etc.\n"
    // uncomment to document debug logging:
    //    + "  -log <function_name>\n"
    + "                  Various levels of logging (-nocolor for monochrome)\n"
    + "  -V, --version   Print GNU-style version info\n"
    + "  -h, --help      Print help message\n"
    + "  -seed <n>       Seed random number generator\n"
    + "\n"
    + "REFERENCES\n"
    + "\n"
    + "The reconstruction method uses phylogenetic transducers, as described in:\n"
    + "  Westesson, Lunter, Paten & Holmes (2012). Accurate Reconstruction of\n"
    + "  Insertion-Deletion Histories by Statistical Phylogenetics.\n"
    + "  PLoS One, DOI: 10.1371/journal.pone.0034572\n"
    + "  http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0034572\n"
    + "\n"
    + "A longer, tutorial-style introduction is available here:\n"
    + "  Westesson, Lunter, Paten & Holmes (2012).\n"
    + "  Phylogenetic Automata, Pruning, and Multiple Alignment.\n"
    + "  http://arxiv.org/abs/1103.4347\n"
    + "\n"
    + "Model-fitting uses the following phylogenetic EM algorithm:\n"
    + "  Holmes & Rubin (2002). An Expectation Maximization Algorithm\n"
    + "  for Training Hidden Substitution Models.\n"
    + "  Journal of Molecular Biology, 317(5).\n"
    + "\n";
}

int main (int argc, char** argv) {

  ProgUsage usage (argc, argv);
  Reconstructor recon;

  const string command = usage.getCommand();
  deque<string>& argvec (usage.argvec);

  if (command == "reconstruct" || command == "recon" || command == "r") {

    recon.reconstructRoot = true;
    recon.accumulateSubstCounts = false;
    recon.accumulateIndelCounts = false;

    usage.implicitSwitches.push_back (string ("-auto"));
    usage.unlimitImplicitSwitches = true;

    while (logger.parseLogArgs (argvec)
	   || recon.parseReconArgs (argvec)
	   || usage.parseUnknown())
      { }

    recon.loadModel();
    recon.loadSeqs();

    recon.reconstructAll();
    recon.writeRecon (cout);
    
  } else if (command == "count" || command == "c") {

    recon.reconstructRoot = false;
    recon.accumulateSubstCounts = true;
    recon.accumulateIndelCounts = true;

    usage.implicitSwitches.push_back (string ("-auto"));
    usage.unlimitImplicitSwitches = true;

    while (logger.parseLogArgs (argvec)
	   || recon.parseFitArgs (argvec)
	   || usage.parseUnknown())
      { }

    recon.loadModel();
    recon.loadSeqs();
    recon.loadRecon();
    recon.loadCounts();

    recon.countAll();
    recon.writeCounts (cout);

  } else if (command == "sum" || command == "s") {

    recon.useLaplacePseudocounts = false;
    
    usage.implicitSwitches.push_back (string ("-counts"));
    usage.unlimitImplicitSwitches = true;
    
    while (logger.parseLogArgs (argvec)
	   || recon.parseSumArgs (argvec)
	   || usage.parseUnknown())
      { }

    recon.loadCounts();
    recon.writeCounts (cout);

  } else if (command == "fit" || command == "f") {

    recon.reconstructRoot = false;
    recon.accumulateSubstCounts = true;
    recon.accumulateIndelCounts = true;

    usage.implicitSwitches.push_back (string ("-auto"));
    usage.unlimitImplicitSwitches = true;
    
    while (logger.parseLogArgs (argvec)
	   || recon.parseFitArgs (argvec)
	   || usage.parseUnknown())
      { }

    recon.loadModel();
    recon.loadSeqs();
    recon.loadRecon();
    recon.loadCounts();
    recon.fit();
    recon.writeModel (cout);
    
  } else
    return usage.parseUnknownCommand (command, HISTORIAN_VERSION);

  return EXIT_SUCCESS;
}
