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
  : OptParser (argc, argv, HISTORIAN_PROGNAME, "{recon,count,fit,mcmc,generate,help,version} [options]")
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
    + "Simulation:\n"
    + "  " + prog + " generate [-model model.json] [-rootlen N] tree.nh >sim.stk\n"
    + "\n"
    + "Commands can be abbreviated to single letters, like so:\n"
    + "  " + prog + " r seqs.fa >reconstruction.stk\n"
    + "  " + prog + " c seqs.fa >counts.json\n"
    + "  " + prog + " f -counts counts.json >model.json\n"
    + "  " + prog + " g tree.nh >sim.stk\n"
    + "(etc.)\n"
    + "\n"
    + "If a command is omitted, 'reconstruct' is assumed, e.g.\n"
    + "  " + prog + " data/gp120.fa -v2\n"
    + "\n"
    + "OPTIONS\n"
    + "\n"
    + "Model specification options\n"
    + "~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
    + "  -model <file>   Load substitution & indel model from file (JSON)\n"
    + "  -preset <name>  Select preset model by name\n"
    + "                   (jc, jcrna, dayhoff, jtt, wag, lg, ECMrest, ECMunrest)\n"
    + "\n"
    + "  -normalize      Normalize expected substitution rate\n"
    + "  -insrate <R>, -delrate <R>, -insextprob <P>, -delextprob <P>\n"
    + "                  Override indel parameters\n"
    + "  -inslen <L>, -dellen <L>\n"
    + "                  Alternate way of setting -insextprob & -delextprob\n"
    + "  -gaprate <R>, -gapextprob <P>, -gaplen <L>\n"
    + "                  Shorthand to set both insertion & deletion params\n"
    + "  -subscale <N>, -indelscale <N>, -scale <N>\n"
    + "                  Scale substitution rates, indel rates, or both\n"
    + "\n"
    + "  -gamma <N>      Add N discretized-gamma rate categories\n"
    + "  -shape <S>      Specify shape parameter for gamma distribution\n"
    + "\n"
    + "  -savemodel <f>  Save model to file, prior to any model-fitting\n"
    + "\n"
    + "Reconstruction file I/O options\n"
    + "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
    + "  -auto <file>    Auto-detect file format and guess its purpose\n"
    + "  -seqs <file>    Specify unaligned sequence file (FASTA)\n"
    + "  -guide <file>   Specify guide alignment file (gapped FASTA)\n"
    + "  -tree <file>    Specify phylogeny file (New Hampshire)\n"
    + "  -nexus <file>, -stockholm <file>\n"
    + "                  Specify phylogeny & guide alignment together\n"
    + "\n"
    + "  -saveguide <f>  Save guide alignment to file\n"
    + "                   (guide tree too, if output format allows)\n"
    + "  -output (nexus|fasta|stockholm|json)\n"
    + "                  Specify output format (default is Stockholm)\n"
    + "  -noancs         Do not display ancestral sequences\n"
    + "\n"
    + "  -codon          Interpret sequences as spliced protein-coding DNA/RNA\n"
    + "\n"
    + "Reconstruction algorithm options\n"
    + "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
    + "The reconstruction algorithm iterates through the guide tree in postorder,\n"
    + "aligning each sibling pair and reconstructing a profile of their parent.\n"
    + "The dynamic programming is constrained to a band around a guide alignment.\n"
    + "\n"
    + "  -band <n>       Size of band around guide alignment (default " + to_string(DefaultMaxDistanceFromGuide) + ")\n"
    + "  -noband         Unlimit band, removing dependence on guide alignment\n"
    + "\n"
    + "The reconstructed parent profile is a weighted finite-state transducer\n"
    + "sampled from the posterior distribution implied by the children. The\n"
    + "posterior probability threshold for inclusion in the parent profile and\n"
    + "max number of states in the parent profile can both be specified to trade\n"
    + "sensitivity vs performance.\n"
    + "\n"
    + "  -profminpost <P>, -profsamples <N>\n"
    + "                  Specify minimum posterior prob. (P) for retaining DP states\n"
    + "                   in profile, or sample N traces randomly (default is -profsamples " + to_string(DefaultProfileSamples) + "\n"
    + "  -profmaxstates <S>, -profmaxmem <M>\n"
    + "                  Limit profile to at most S states, or to use at most M% of\n"
    + "                   memory for DP matrix (default is -profmaxmem " + to_string(DefaultMaxDPMemoryFraction) + ")\n"
    //    + "  -profminlen <L>, -profmaxlen <L>\n"
    //    + "                  Constrain permissible range of ancestral sequence lengths\n"
    //    + "                   (use with care; extreme/unreachable values may cause program to hang!)\n"
    + "\n"
    + "Following alignment, ancestral sequence reconstruction can be performed.\n"
    + "\n"
    + "  -ancseq         Predict ancestral sequences (default is to leave them as " + Alignment::wildcardChar + "'s)\n"
    + "  -ancprob        Report posterior probabilities for ancestral residues\n"
    + "\n"
    + "For additional accuracy in historical reconstruction, the alignment can be\n"
    + "iteratively refined, or MCMC-sampled. By default, refinement and MCMC are\n"
    + "both disabled. (MCMC currently requires an ultrametric tree.)\n"
    + "\n"
    + "  -norefine, -refine"
    + "                  Disable/enable iterative refinement after initial reconstruction\n"
    + "\n"
    + "  -mcmc           Run MCMC sampler after reconstruction\n"
    + "  -samples <N>    Number of MCMC iterations per sequence (default " + to_string(DefaultMCMCSamplesPerSeq) + ")\n"
    + "  -trace <file>   Specify MCMC trace filename\n"
    + "  -fixtree        Fix tree during MCMC (sample alignment only)\n"
    + "  -fixalign       Fix alignment during MCMC (sample tree only)\n"
    //    + "  -fixguide       Fix guide alignment during MCMC\n"
    + "\n"
    + "Guide alignment & tree estimation options\n"
    + "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
    + "The guide aligner builds a maximal spanning tree of pairwise alignments.\n"
    + "It can be accelerated in two ways. The first is by using a sparse random\n"
    + "graph instead of a fully connected all-vs-all pairwise comparison.\n"
    + "\n"
    + "  -rndspan        Use a sparse random spanning graph (default)\n"
    + "  -allspan        Use a dense random spanning graph, i.e. all-vs-all pairs\n"
    + "\n"
    + "The second way to optimize construction of the guide alignment is by\n"
    + "confining the pairwise DP matrix to cells around a subset of diagonals\n"
    + "that contain above a threshold number of k-mer matches. To turn on the\n"
    + "former optimization, use -rndspan; the latter is turned on by default for\n"
    + "sequences whose full DP matrix would not otherwise fit in memory (the\n"
    + "memory threshold can be set with -kmatchmb). It can be disabled with\n"
    + "-kmatchoff, or enabled (for a particular k-mer threshold) with -kmatchn.\n"
    + "\n"
    + "  -kmatchn <n>    Threshold# of kmer matches to seed a diagonal\n"
    + "                   (default sets this as low as available memory will allow)\n"
    + "  -kmatch <k>     Length of kmers for pre-filtering heuristic (default " + to_string(DEFAULT_KMER_LENGTH) + ")\n"
    + "  -kmatchband <n> Size of DP band around kmer-matching diagonals (default " + to_string(DEFAULT_BAND_SIZE) + ")\n"
    + "  -kmatchmb <M>   Set kmer threshold to use M megabytes of memory\n"
    + "  -kmatchmax      Set kmer threshold to use all available memory\n"
    + "  -kmatchoff      No kmer threshold, do full DP\n"
    + "\n"
    + "Following construction of the guide alignment, a tree is estimated using a\n"
    + "distance matrix method. By default this is UPGMA.\n"
    + "\n"
    + "  -upgma          Use UPGMA to estimate tree (default for MCMC)\n"
    + "  -nj             Use neighbor-joining, not UPGMA, to estimate tree\n"
    + "  -jc             Use Jukes-Cantor-like estimates for distance matrix\n"
    + "\n"
    + "Some common settings (the default is somewhere in between these extremes):\n"
    + "\n"
    + "  -careful        Run in careful mode. Shorthand for the following:\n"
    + "                   " + join(Reconstructor::carefulAliasArgs," ") + "\n"
    + "\n"
    + "  -fast           Run in fast mode. Shorthand for the following:\n"
    + "                   " + join(Reconstructor::fastAliasArgs," ") + "\n"
    + "\n"
    + "Model-fitting and event-counting options\n"
    + "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
    + "In reconstruction mode, any supplied alignment will be interpreted as a hint,\n"
    + "i.e. a guide alignment. In MCMC, counting, or model-fitting mode, any alignment\n"
    + "that contains a full ancestral sequence reconstruction will be interpreted as a\n"
    + "reconstruction. To force the alignment to be interpreted as a reconstruction,\n"
    + "precede it with -recon, -nexusrecon or -stockrecon (depending on the format).\n"
    + "\n"
    + "  -recon <file>, -nexusrecon <file>, -stockrecon <file>\n"
    + "                  Use precomputed reconstruction (FASTA/NEXUS/Stockholm)\n"
    + "  -mininc <n>     EM convergence threshold as relative log-likelihood increase\n"
    + "                   (default is " + TOSTRING(DefaultMinEMImprovement) + ")\n"
    + "  -maxiter <n>    Max number of EM iterations (default " + to_string(DefaultMaxEMIterations) + ")\n"
    + "  -nolaplace      Do not add Laplace +1 pseudocounts during model-fitting\n"
    + "  -fixsubrates    Do not estimate substitution rates or initial composition\n"
    + "  -fixgaprates    Do not estimate indel rates or length distributions\n"
    + "\n"
    + "General options\n"
    + "~~~~~~~~~~~~~~~\n"
    + "  -verbose, -vv, -vvv, -v4, -v5, etc.\n"
    // uncomment to document debug logging:
    //    + "  -log <function_name>\n"
    + "                  Various levels of logging (-nocolor for monochrome)\n"
    + "  -V, --version   Print GNU-style version info\n"
    + "  -h, --help      Print help message\n"
    + "  -seed <n>       Seed random number generator (" + DPMatrix::random_engine_name() + "; default seed " + to_string(DPMatrix::random_engine::default_seed) + ")\n"
    + "\n"
    + "REFERENCES\n"
    + "\n"
    + "The primary reference for this software is the following:\n"
    + "  Holmes (2017). Historian: Accurate Reconstruction of\n"
    + "  Ancestral Sequences and Evolutionary Rates.\n"
    + "  Bioinformatics, DOI: 10.1093/bioinformatics/btw791\n"
    + "  https://academic.oup.com/bioinformatics/article/33/8/1227/2926463\n"
    + "\n"
    + "The reconstruction algorithm uses phylogenetic transducers, as described in:\n"
    + "  Westesson, Lunter, Paten & Holmes (2012). Accurate Reconstruction of\n"
    + "  Insertion-Deletion Histories by Statistical Phylogenetics.\n"
    + "  PLoS One, DOI: 10.1371/journal.pone.0034572\n"
    + "  http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0034572\n"
    + "\n"
    + "A longer, tutorial-style introduction to transducers is available here:\n"
    + "  Westesson, Lunter, Paten & Holmes (2012).\n"
    + "  Phylogenetic Automata, Pruning, and Multiple Alignment.\n"
    + "  http://arxiv.org/abs/1103.4347\n"
    + "\n"
    + "Model-fitting uses the following phylogenetic EM algorithm:\n"
    + "  Holmes & Rubin (2002). An Expectation Maximization Algorithm\n"
    + "  for Training Hidden Substitution Models.\n"
    + "  Journal of Molecular Biology, 317(5).\n"
    + "\n"
    + "The MCMC kernels for co-sampling alignments and trees are described in:\n"
    + "  Holmes & Bruno (2001). Evolutionary HMMs: A Bayesian Approach to\n"
    + "  Multiple Alignment. Bioinformatics, 17(9).\n"
    + "  Redelings & Suchard (2005). Joint Bayesian Estimation of Alignment\n"
    + "  and Phylogeny. Systematic Biology, 54(3).\n"
    + "\n";
}

int main (int argc, char** argv) {

  ProgUsage usage (argc, argv);
  Reconstructor recon;

  const string command = usage.getCommand();
  deque<string>& argvec (usage.argvec);

  auto reconstruct = [&]() -> void
    {
      recon.reconstructRoot = true;
      recon.accumulateSubstCounts = false;
      recon.accumulateIndelCounts = false;

      usage.implicitSwitches.push_back (string ("-auto"));
      usage.unlimitImplicitSwitches = true;

      while (logger.parseLogArgs (argvec)
	     || recon.parseReconArgs (argvec)
	     || recon.parseModelArgs (argvec)
	     || recon.parseProfileArgs (argvec, false)
	     || recon.parseSamplerArgs (argvec)
	     || recon.parseAncSeqArgs (argvec)
	     || usage.parseUnknown())
	{ }

      recon.loadModel();
      recon.loadSeqs();

      recon.reconstructAll();
      recon.sampleAll();
      recon.predictAllAncestors();
      
      recon.writeRecon (cout);
    };

  if (command == "reconstruct" || command == "recon" || command == "r") {
    reconstruct();

  } else if (command == "generate" || command == "g") {

    usage.implicitSwitches.push_back (string ("-tree"));
    usage.unlimitImplicitSwitches = true;

    while (logger.parseLogArgs (argvec)
	   || recon.parseSimulatorArgs (argvec)
	   || recon.parseModelArgs (argvec)
	   || usage.parseUnknown())
      { }

    recon.simulate();
    
  } else if (command == "mcmc" || command == "m") {

    recon.reconstructRoot = true;
    recon.runMCMC = true;
    recon.useUPGMA = true;
    recon.accumulateSubstCounts = false;
    recon.accumulateIndelCounts = false;

    usage.implicitSwitches.push_back (string ("-auto"));
    usage.unlimitImplicitSwitches = true;

    while (logger.parseLogArgs (argvec)
	   || recon.parsePremadeArgs (argvec)
	   || recon.parseModelArgs (argvec)
	   || recon.parseProfileArgs (argvec, true)
	   || recon.parseSamplerArgs (argvec)
	   || recon.parseAncSeqArgs (argvec)
	   || usage.parseUnknown())
      { }

    recon.loadModel();
    recon.loadSeqs();
    recon.loadRecon();

    recon.sampleAll();
    recon.writeRecon (cout);

  } else if (command == "count" || command == "c") {

    recon.reconstructRoot = false;
    recon.accumulateSubstCounts = true;
    recon.accumulateIndelCounts = true;

    usage.implicitSwitches.push_back (string ("-auto"));
    usage.unlimitImplicitSwitches = true;

    while (logger.parseLogArgs (argvec)
	   || recon.parsePremadeArgs (argvec)
	   || recon.parseModelArgs (argvec)
	   || recon.parseProfileArgs (argvec, true)
	   || recon.parseCountArgs (argvec)
	   || recon.parseSumArgs (argvec)
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
	   || recon.parsePremadeArgs (argvec)
	   || recon.parseModelArgs (argvec)
	   || recon.parseProfileArgs (argvec, true)
	   || recon.parseCountArgs (argvec)
	   || recon.parseSumArgs (argvec)
	   || recon.parseFitArgs (argvec)
	   || usage.parseUnknown())
      { }

    recon.loadModel();
    recon.loadSeqs();
    recon.loadRecon();
    recon.loadCounts();
    recon.fit();
    recon.writeModel (cout);
    
  } else if (!usage.parseUnknownCommand (command, HISTORIAN_VERSION, false)) {

    // default: reconstruct
    argvec.push_front (command);
    reconstruct();
  }

  return EXIT_SUCCESS;
}
