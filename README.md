# Indel Historian

Indel Historian is a multiple aligner that aims at providing accurate historical reconstructions of the evolution of a set of DNA or protein sequences. This stands in contrast to most multiple alignment tools, which instead try to provide protein alignments that correctly identify regions of 3D structural homology. If you are trying to predict the structure of a protein, you should probably use one of those other tools. If you care about the evolutionary history of your sequences, use Indel Historian.

Most multiple alignment tools (think: Clustal, Muscle, ProbCons) are optimized for homology-based structure prediction, and tested on structural alignment benchmarks (e.g. BAliBase, Oxbench, Prefab, Sabmark). That's a good empirical approach as far as it goes, because structurally-informed protein alignments make a good "gold standard" for benchmarking alignment tools. Often, these tools have scoring schemes that are optimized for reproducing common signatures of protein selection, such as reduced indel rates in hydrophobic regions. However, optimizing for structure has the unfortunate drawback of introducing biases into the estimates of indel (and possibly substitution) rates (as a rule of thumb, all methods tend to underestimate mutation rates, but with standard multiple aligners the biases can be unpredictable and can vary widely at different indel rates). Consequently, this approach yields a less-than-accurate picture of evolutionary history.

By contrast, Indel Historian uses an explicit evolutionary model of indel and substitution events, derived rigorously from statistical phylogenetics using finite-state transducers as evolutionary operators. In simulation tests (using the third-party evolution simulator [indel-Seq-Gen](https://www.ncbi.nlm.nih.gov/pubmed/17158778)), it introduces significantly fewer biases than other tools. It also performs pretty well on structural alignment benchmarks, though not as well as tools like Muscle and ProbCons that are optimized for that.

The basic method and initial benchmarks of the approach are described in ([_Accurate Reconstruction of Insertion-Deletion Histories by Statistical Phylogenetics_](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0034572); Westesson, Lunter, Paten and Holmes, PLoS One, 2012). The implementation described in that paper is a program called ProtPal. Indel Historian is a clean reimplementation of ProtPal that also runs a lot faster, is more user-friendly, and has significantly more features.

The closest method to Indel Historian is PRANK ([_Phylogeny-Aware Gap Placement Prevents Errors in Sequence Alignment and Evolutionary Analysis_](https://www.ncbi.nlm.nih.gov/pubmed/18566285); Löytynoja and Goldman, Science, 2008). Indel Historian is slightly more accurate than PRANK (on both simulation and structural alignment benchmarks), has more features, and runs faster.

## Installation

To build from source, type `make`. This will create a binary file `bin/historian`.

At present, Indel Historian requires the following build environment to compile:

* Apple LLVM version 7.3.0 or later (clang-703.0.31), or gcc version 4.8.3 or later
 * Boost C++ library version 1.62.0 or later
 * Gnu Scientific Library (GSL) version 2.2.1 or later

If you want to run the tests, type `make test`. (The tests are oriented to a Mac OSX build environment; numerical precision errors may cause slight differences in output on different machines, which may lead to some tests failing.)

Pre-compiled binaries are also available from the GitHub repository [release page](https://github.com/ihh/indelhistorian/releases).

## Examples

### Basic reconstruction

The simplest way to use Indel Historian is just to point it at a FASTA file. It will then estimate a guide alignment, estimate a tree from that (using neighbor-joining), and perform a full ancestral reconstruction.

For example, using a test file of [HIV GP120 sequences](https://github.com/ihh/indelhistorian/blob/master/data/gp120.fa) that is included in the repository:

	historian data/gp120.fa

This will generally be pretty fast, but you can make it faster (at a slight cost in accuracy) using the `-fast` option:

	historian data/gp120.fa -fast

The `-fast` option is an alias for several reconstruction options, as described in the [help message](#HelpText).

#### Commands

Indel Historian is one of these toolbox programs where the first argument can be a command, specifying what action is to be performed. If you omit this first command-argument, Indel Historian assumes you want to reconstruct something. You can make this explicit as follows:

	historian reconstruct data/gp120.fa

You can also abbreviate `reconstruct` to `recon` or just `r`:

	historian r data/gp120.fa

All commands can be abbreviated to single letters.
Of course, if you are doing reconstruction (as noted above), you can omit the command entirely, but the examples [further below](#ModelFitting) will use different commands.

#### Logging

By default, the `historian` program runs in classic Unix mode, i.e. printing nothing except the output. This can be a bit boring for long jobs, so run with `-v` for verbose colorized logging output, or `-vv` for more verbose, or `-vvv` (or equivalently `-v3`) for even more logging, and so on. Probably `-v2` is about the right balance between showing some feedback on progress and not drowning you in psychedelic rainbow log messages.

	historian -v2 data/gp120.fa

This produces output somewhat like this (but more colorful - GitHub-flavored Markdown won't show the colors):

	Using default amino acid model
	Loading sequences from data/gp120.fa
	Building guide alignment (data/gp120.fa)
	Estimating initial tree by neighbor-joining (data/gp120.fa)
	Starting reconstruction on 19-node tree (data/gp120.fa)
	Aligning ENV_HV2BE/24-510 (489 states, 488 transitions) and ENV_HV2D1/24-501 (480 states, 479 transitions)
	Aligning ENV_HV2G1/23-502 (482 states, 481 transitions) and (ENV_HV2BE/24-510:0.106114,ENV_HV2D1/24-501:0.0784677) (689 states, 808 transitions)
...and so on.

#### File formats

Indel Historian speaks a variety of input and output formats. By default, it outputs alignments in [Stockholm format](https://en.wikipedia.org/wiki/Stockholm_format), which allows easy extraction of the alignment while also affording space for metadata like trees. If you prefer your alignments in another output format, for example [Nexus](https://en.wikipedia.org/wiki/Nexus_file) or [FASTA](https://en.wikipedia.org/wiki/FASTA_format), use `-output nexus` or `-output fasta`.

Input sequence formats will usually be auto-detected, but this behavior can be overridden to stipulate particular file formats via the [command-line arguments](#HelpText).

#### Fine-tuning the reconstruction

As briefly alluded to above, Indel Historian does several performance-optimizing steps _en route_ to a reconstruction. First, it builds a quick-guess multiple alignment by a greedy maximal-spanning-tree type approach; this can optionally be accelerated by a k-mer match step (confining the alignment to diagonals of the dynamic programming matrix that pass a minimum threshold of k-mer matches) and by using a sparse [random spanning forest](https://www.ncbi.nlm.nih.gov/pubmed/19478997) instead of a dense all-vs-all comparison. Second, it uses this alignment to build a guide tree by neighbor-joining. Third, it builds a progressive reconstruction that includes suboptimal alignments in something like a [partial-order graph](https://www.ncbi.nlm.nih.gov/pubmed/11934745). And fourth, it optionally does iterative refinement to optimize the reconstruction. The latter two steps (reconstruction and refinement) can be constrained to stay near the guide alignment for performance reasons.

The default settings attempt to navigate this maze of options for you, mostly using the higher-accuracy options until memory becomes a limiting factor and then switching to the more approximate options. However, as a power user, you may want to take control of these options. Command-line arguments allow you to supply guide alignments and/or guide trees, and change the parameters or behavior of the standard workflow.

For example, if you want to supply sequences in FASTA format ([gp120.fa](https://github.com/ihh/indelhistorian/blob/master/data/gp120.fa)) and a guide tree in New Hampshire format ([gp120.tree.nh](https://github.com/ihh/indelhistorian/blob/master/data/gp120.tree.nh))

	historian -seqs data/gp120.fa -tree gp120.tree.nh

Alternatively, if your sequences are in Nexus or Stockholm format, you can encode the tree together with your sequences, using the appropriate syntax for encoding New Hampshire-format trees in those formats (`#=GF NH` for Stockholm).

If you want to estimate the tree and use UPGMA instead of neighbor-joining, so as to enforce an ultrametric tree, use the `-upgma` option, e.g.

	historian -upgma data/gp120.fa

To save the guide alignment, use the `-saveguide` option:

	historian data/gp120.fa -saveguide gp120.guide.fa

If you already have your sequences aligned, and you want to use [this alignment](https://github.com/ihh/indelhistorian/blob/master/data/gp120.guide.fa) as the guide alignment, you only need to supply that alignment (not the sequences):

	historian -guide gp120.guide.fa

The guide alignment is, by default, just treated as a "hint". Indel Historian will do dynamic programming in a "band" around the guide alignment, sliding gaps back and forth up to a maximum distance specified by the `-band` argument. If, instead, you want to use the guide alignment as a strict constraint, and find the best reconstruction that is exactly consistent with the guide, then set the band to zero:

	historian -guide gp120.guide.fa -band 0

For some alignments, it may be the case that there is no reconstruction under Indel Historian's model that is exactly consistent with the guide (due to ordering of gaps), in which case you might want to relax the band to 1:

	historian -guide gp120.guide.fa -band 1

These arguments are all listed in the help text, available via the `-h` option and copied [below](#HelpText).

## Model-fitting

Indel Historian's underlying model is a simple one: there is a substitution rate matrix, an insertion rate, a deletion rate, and insertion/deletion extension probabilities. These are all specified in a JSON file format, an example of which can be found [here](https://github.com/ihh/indelhistorian/blob/master/data/testcount.jukescantor.json).

The default model `lg` is an amino acid substitution matrix estimated by [Le and Gascuel (2008)](https://www.ncbi.nlm.nih.gov/pubmed/18367465) using [XRate](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0036898) on a dataset of [Pfam](http://pfam.xfam.org/) alignments, with indel rates and probabilities that were also estimated from Pfam. However, the `historian` program allows you to load a model from a file using the `-model` option, or to use one of the preset models using the `-preset` option. Alternatively, the model parameters can be estimated directly from sequence data using the built-in EM algorithm that is the same algorithm used by XRate (as described in [this paper](http://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-7-428)).

To estimate rates from data, use the `fit` command. Model-fitting takes a little longer than reconstruction, since the EM algorithm typically takes a few iterations to converge, so you might want to turn on some logging. For example:

	historian fit data/gp120.fa -fast -vv >gp120.model.json

You can then load this model using the `-model` option, and use it to reconstruct another sequence history. For example, using it to reconstruct the ancestors of the Cas9 bridge helix domain (Pfam family PF16593), included as file [data/PF16593.fa](https://github.com/ihh/indelhistorian/blob/master/data/PF16593.fa) in the repository:

	historian reconstruct -model gp120.model.json data/PF16593.fa

Often, a single protein family will not include enough information to reliably fit an amino acid rate matrix - the model will be over-trained. You can aggregate multiple datasets together into a larger training set simply by listing the files successively on the command line:

	historian fit -fast data/gp120.fa data/PF16593.fa -v3 >aggregated.model.json

If you only want to estimate the indel rates and not the substitution matrix, then you can use the `-fixsubrates` option to hold the substitution rates constant:

	historian fit data/gp120.fa -fixsubrates >gp120.model.json

Conversely, you can use `-fixgaprates` to hold the indel rates (and indel extension parameters) constant, while estimating substitution rates. Other aspects of the model-fitting algorithm (for example, the use of [Laplace pseudocounts](https://en.wikipedia.org/wiki/Additive_smoothing), or the EM convergence criteria) can be set via the [command-line options](#HelpText).

## Event-counting

The model-fitting algorithm described above is an [Expectation-Maximization algorithm](https://en.wikipedia.org/wiki/Expectation%E2%80%93maximization_algorithm) that iteratively estimates the _a posteriori_ expected number of indel and substitution events of each possible type, along with the amount of time spent in each state, and then uses these expected counts and times to update the mutation rates.

For some applications, these event counts and wait times - the so-called "sufficient statistics" for EM - may be more useful (or more stable, or reliable) than the rate parameters themselves.

If you are just interested in these counts, you can get at them using the `count` command. For example:

	historian count data/gp120.fa >gp120.counts.json

One application of the counts is as [pseudocounts](https://en.wikipedia.org/wiki/Pseudocount) representing hyperparameters of Dirichlet priors. JSON counts files that were computed using the `count` command can be passed as pseudocounts to the `fit` command via the `-counts` option:

	historian fit data/PF16593.fa -counts gp120.counts.json -nolaplace >PF16593.model.json

The `-nolaplace` option here indicates that we don't want to add the usual +1 Laplace pseudocounts in this case.

You can also use the `count` command to estimate counts for a bunch of alignments in parallel, then combine the counts using the `sum` command, and finally estimate the maximum-likelihood parameters for these combined counts using `fit`. This amounts to doing a single iteration of the EM algorithm. Probably the only reason you'd want to do this would be if you were implementing some kind of roll-your-own [MapReduce](https://en.wikipedia.org/wiki/MapReduce)-style approach to splitting up the EM algorithm on a large dataset. Which is probably unlikely, but here's roughly what it'd look like:

	historian count file1.fa >file1.counts.json
	historian count file2.fa >file2.counts.json
	...
	historian count fileN.fa >fileN.counts.json
	historian sum file1.counts.json file2.counts.json ... fileN.counts.json >summed.counts.json
	historian fit -nolaplace -counts summed.counts.json >updated.model.json

## MCMC

Indel Historian includes an experimental MCMC implementation for co-sampling trees and alignments. Currently, this implementation only works for ultrametric trees. It is available via the `mcmc` command.

## Method
At its core, Indel Historian uses the phylogenetic transducer method.
See [Westesson et al, 2012](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0034572) for an evaluation and brief description of the method, or [this arXiv report](http://arxiv.org/abs/1103.4347) for a tutorial introduction.

Very briefly, the idea of this method is as follows. The main recursion of Felsenstein's [pruning algorithm](https://en.wikipedia.org/wiki/Felsenstein%27s_tree-pruning_algorithm) for calculating the likelihood of a multiple alignment column can be summarized, in matrix form, as **Fn=(Bl Fl).(Br Fr)** where **n**, **l** and **r** are the node and its two children, **Bn** is the branch substitution matrix on the branch leading to node **n**, **(A B)** denotes the matrix product and **A.B** the pointwise (Hadamard) product, with each **Fn** denoting an ancestral sequence profile. If for our matrix representation we use [weighted finite-state transducers](https://en.wikipedia.org/wiki/Finite-state_transducer), with **(A B)** denoting the operation of transducer composition and **A.B** the operation of transducer intersection, then Felsenstein's algorithm yields an instance of [Sankoff's algorithm](http://epubs.siam.org/doi/abs/10.1137/0145048) for multiple sequence alignment, and **Fn** is a state machine. We constrain the algorithm to be practical by retaining only high-probability states of **Fn** at each stage. The branch transducers **Bn** are derived using a simple approximation that indel events on a single branch never overlap.

## Release Notes, v1.0 (10/9/2016)

* The tests work for the following development environment. Other builds may give subtly different results due to rounding errors.
 * Apple LLVM version 7.3.0 (clang-703.0.31)
 * Boost 1.62.0
 * GSL 2.2.1
* At present, IndelHistorian requires the Boost C++ library due to inconsistencies in the way different C++ compilers and standard libraries implement regular expressions.

<a name="HelpText"></a>
## Command-line help text

The following is the message that appears when you type `historian help`:

<pre><code>
Usage: historian {recon[struct],count,fit,mcmc,generate,help,version} [options]

EXAMPLES

Reconstruction:
  historian recon seqs.fa [-tree tree.nh] -output fasta &gt;reconstruction.fa
  historian recon -guide guide.fa [-tree tree.nh] &gt;reconstruction.stk
  historian recon guide.stk &gt;reconstruction.stk
  historian recon data.nex -output nexus &gt;reconstruction.nex

Event counting:
  historian count seqs.fa [-tree tree.nh] [-model model.json] &gt;counts.json
  historian count -guide guide.fa [-tree tree.nh] &gt;counts.json
  historian count -recon reconstruction.fa -tree tree.nh &gt;counts.json

Model fitting:
  historian fit seqs.fa &gt;newmodel.json
  historian fit -counts counts.json &gt;newmodel.json

Simulation:
  historian generate [-model model.json] [-rootlen N] tree.nh &gt;sim.stk

Commands can be abbreviated to single letters, like so:
  historian r seqs.fa &gt;reconstruction.stk
  historian c seqs.fa &gt;counts.json
  historian f -counts counts.json &gt;model.json
  historian g tree.nh &gt;sim.stk
(etc.)

If a command is omitted, 'reconstruct' is assumed.

OPTIONS

Reconstruction file I/O options
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  -auto &lt;file&gt;    Auto-detect file format and guess its purpose
  -model &lt;file&gt;   Specify substitution & indel model file (JSON)
  -preset &lt;name&gt;  Select preset model by name
                   (gtr, dayhoff, jtt, wag, lg, ECMrest, ECMunrest)
  -seqs &lt;file&gt;    Specify unaligned sequence file (FASTA)
  -guide &lt;file&gt;   Specify guide alignment file (gapped FASTA)
  -tree &lt;file&gt;    Specify phylogeny file (New Hampshire)
  -nexus &lt;file&gt;, -stockholm &lt;file&gt;
                  Specify phylogeny & guide alignment together

  -saveguide &lt;f&gt;  Save guide alignment to file
                   (guide tree too, if output format allows)
  -output (nexus|fasta|stockholm)
                  Specify output format (default is Stockholm)
  -noancs         Do not display ancestral sequences

  -codon          Interpret sequences as spliced protein-coding DNA/RNA

Reconstruction algorithm options
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The reconstruction algorithm iterates through the guide tree in postorder,
aligning each sibling pair and reconstructing a profile of their parent.
The dynamic programming is constrained to a band around a guide alignment.

  -band &lt;n&gt;       Size of band around guide alignment (default 40)
  -noband         Unlimit band, removing dependence on guide alignment

The reconstructed parent profile is a weighted finite-state transducer
sampled from the posterior distribution implied by the children. The
posterior probability threshold for inclusion in the parent profile and
max number of states in the parent profile can both be tweaked to trade off
sensitivity vs performance.

  -profminpost &lt;P&gt;, -profsamples &lt;N&gt;
                  Specify minimum posterior prob. (P) for retaining DP states
                   in profile (default .01), or sample N traces randomly
  -profmaxstates &lt;S&gt;
                  Limit profile to at most S states
                   (default uses at most 5% of memory for DP matrix)

Following alignment, ancestral sequence reconstruction can be performed.

  -ancseq         Predict ancestral sequences (default is to leave them as *'s)
  -ancprob        Report posterior probabilities for ancestral residues

For additional accuracy in historical reconstruction, the alignment can be
iteratively refined, or MCMC-sampled. By default, refinement is enabled and
MCMC is disabled. (MCMC currently requires an ultrametric tree.)

  -norefine       Disable iterative refinement after initial reconstruction

  -mcmc           Run MCMC sampler after reconstruction
  -samples &lt;N&gt;    Number of MCMC iterations per sequence (default 100)
  -trace &lt;file&gt;   Specify MCMC trace filename

Guide alignment & tree estimation options
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The guide aligner builds a maximal spanning tree of pairwise alignments.
It can be accelerated in two ways. The first is by using a sparse random
forest instead of a fully connected all-vs-all pairwise comparison.

  -rndspan        Use a sparse random spanning graph, not all-vs-all pairs

The second way to optimize construction of the guide alignment is by
confining the pairwise DP matrix to cells around a subset of diagonals
that contain above a threshold number of k-mer matches. To turn on the
former optimization, use -rndspan; the latter is turned on by default for
sequences whose full DP matrix would not otherwise fit in memory (the
memory threshold can be set with -kmatchmb). It can be disabled with
-kmatchoff, or enabled (for a particular k-mer threshold) with -kmatchn.

  -kmatchn &lt;n&gt;    Threshold# of kmer matches to seed a diagonal
                   (default sets this as low as available memory will allow)
  -kmatch &lt;k&gt;     Length of kmers for pre-filtering heuristic (default 6)
  -kmatchband &lt;n&gt; Size of DP band around kmer-matching diagonals (default 64)
  -kmatchmb &lt;M&gt;   Set kmer threshold to use M megabytes of memory
  -kmatchoff      No kmer threshold, do full DP

Following construction of the guide alignment, a tree is estimated using a
distance matrix method. By default this is neighbor-joining.

  -upgma          Use UPGMA to estimate tree (default for MCMC)
  -nj             Use neighbor-joining, not UPGMA, to estimate tree
  -jc             Use Jukes-Cantor-like estimates for distance matrix

If you are confident the guide alignment & tree should be reasonably obvious,
and just want to get on to reconstruction as quickly as possible:

  -fast           Shorthand for the following options:
                   -rndspan -kmatchn 3 -band 10 -profmaxstates 1 -jc -norefine

Model-fitting and event-counting options
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
In reconstruction mode, any supplied alignment will be interpreted as a hint,
i.e. a guide alignment. In MCMC, counting, or model-fitting mode, any alignment
that contains a full ancestral sequence reconstruction will be interpreted as a
reconstruction. To force the alignment to be interpreted as a reconstruction,
precede it with -recon, -nexusrecon or -stockrecon (depending on the format).

  -recon &lt;file&gt;, -nexusrecon &lt;file&gt;, -stockrecon &lt;file&gt;
                  Use precomputed reconstruction (FASTA/NEXUS/Stockholm)
  -mininc &lt;n&gt;     EM convergence threshold as relative log-likelihood increase
                   (default is .001)
  -maxiter &lt;n&gt;    Max number of EM iterations (default 100)
  -nolaplace      Do not add Laplace +1 pseudocounts during model-fitting
  -fixsubrates    Do not estimate substitution rates or initial composition
  -fixgaprates    Do not estimate indel rates or length distributions

General options
~~~~~~~~~~~~~~~
  -verbose, -vv, -vvv, -v4, -v5, etc.
                  Various levels of logging (-nocolor for monochrome)
  -V, --version   Print GNU-style version info
  -h, --help      Print help message
  -seed &lt;n&gt;       Seed random number generator (mt19937; default seed 5489)

REFERENCES

The initial reconstruction uses phylogenetic transducers, as described in:
  Westesson, Lunter, Paten & Holmes (2012). Accurate Reconstruction of
  Insertion-Deletion Histories by Statistical Phylogenetics.
  PLoS One, DOI: 10.1371/journal.pone.0034572
  http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0034572

A longer, tutorial-style introduction is available here:
  Westesson, Lunter, Paten & Holmes (2012).
  Phylogenetic Automata, Pruning, and Multiple Alignment.
  http://arxiv.org/abs/1103.4347

Model-fitting uses the following phylogenetic EM algorithm:
  Holmes & Rubin (2002). An Expectation Maximization Algorithm
  for Training Hidden Substitution Models.
  Journal of Molecular Biology, 317(5).

The MCMC kernels for co-sampling alignments and trees are described in:
  Holmes & Bruno (2001). Evolutionary HMMs: A Bayesian Approach to
  Multiple Alignment. Bioinformatics, 17(9).
  Redelings & Suchard (2005). Joint Bayesian Estimation of Alignment
  and Phylogeny. Systematic Biology, 54(3).

</code></pre>
