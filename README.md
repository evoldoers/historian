# Indel Historian
Reconstruction of evolutionary indel & substitution histories using the phylogenetic transducer method
(see [Westesson et al, 2012](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0034572) for an evaluation and brief description of the method, or [this arXiv report](http://arxiv.org/abs/1103.4347) for a tutorial introduction).

<pre><code>
Usage: historian {recon[struct],count,fit,help,version} [options]

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

Commands can be abbreviated to single letters, like so:
  historian r seqs.fa &gt;reconstruction.stk
  historian c seqs.fa &gt;counts.json
  historian f -counts counts.json &gt;model.json
(etc.)

OPTIONS

Reconstruction file I/O options:
  -auto &lt;file&gt;    Auto-detect file format and guess its purpose
  -model &lt;file&gt;   Specify substitution & indel model (JSON)
  -seqs &lt;file&gt;    Specify unaligned sequences (FASTA)
  -guide &lt;file&gt;   Specify guide alignment (gapped FASTA)
  -tree &lt;file&gt;    Specify phylogeny (New Hampshire)
  -nexus &lt;file&gt;, -stockholm &lt;file&gt;
                  Specify phylogeny & guide alignment together

  -saveguide &lt;f&gt;  Save guide alignment to file
                   (guide tree too, if output format allows)
  -output (nexus|fasta|stockholm)
                  Specify output format (default is Stockholm)

Reconstruction algorithm options:

The reconstruction algorithm iterates through the guide tree in postorder,
aligning each sibling pair and reconstructing a profile of their parent.
The dynamic programming is constrained to a band around a guide alignment,
and the reconstructed parent profile is a weighted finite-state transducer
sampled from the posterior distribution implied by the children. The band
size, posterior probability threshold for inclusion in the parent profile,
and number of states in the parent profile can all be tweaked to trade off
sensitivity vs performance.

  -band &lt;n&gt;       Size of band around guide alignment (default 10)
  -noband         Unlimit band, ignore guide alignment
  -minpost &lt;p&gt;    Posterior prob. threshold for profile states (default .01)
  -states &lt;n&gt;     Limit max number of states per profile

  -ancseq         Predict ancestral sequences (default is to leave them as *'s)

Guide alignment & tree estimation options:

The guide aligner builds a maximal spanning tree of pairwise alignments.
It can be accelerated in two ways: (1) by using a sparse random forest
instead of a fully connected all-vs-all pairwise comparison; and (2) by
confining the pairwise DP matrix to cells around a subset of diagonals
that contain above a threshold number of k-mer matches. To turn on the
former optimization, use -rndspan; the latter is turned on by default for
sequences whose full DP matrix would not otherwise fit in memory (the
memory threshold can be set with -kmatchmb). It can be disabled with
-kmatchoff, or enabled (for a particular k-mer threshold) with -kmatchn.

  -rndspan        Use a sparse random spanning graph, not all-vs-all pairs
  -kmatchn &lt;n&gt;    Threshold# of kmer matches to seed a diagonal
                   (default sets this as low as available memory will allow)
  -kmatch &lt;k&gt;     Length of kmers for pre-filtering heuristic (default 6)
  -kmatchband &lt;n&gt; Size of DP band around kmer-matching diagonals (default 64)
  -kmatchmb &lt;M&gt;   Set kmer threshold to use M megabytes of memory
  -kmatchoff      No kmer threshold, do full DP

  -upgma          Use UPGMA, not neighbor-joining, to estimate tree

Model-fitting and event-counting options:

By default, the reconstruction algorithm will interpret any supplied alignment
as a guide alignment, i.e. a hint, even if it contains a full ancestral sequence
reconstruction. To insist that the alignment be interpreted as a reconstruction,
precede it with -recon, -nexusrecon or -stockrecon (depending on the format).

  -recon &lt;file&gt;, -nexusrecon &lt;file&gt;, -stockrecon &lt;file&gt;
                  Use precomputed reconstruction (FASTA/NEXUS/Stockholm)
  -mininc &lt;n&gt;     EM convergence threshold as relative log-likelihood increase
                   (default is .001)
  -maxiter &lt;n&gt;    Max number of EM iterations (default 100)
  -nolaplace      Do not add Laplace +1 pseudocounts during model-fitting
  -fixsubrates    Do not estimate substitution rates or initial composition
  -fixgaprates    Do not estimate indel rates or length distributions

General options:
  -verbose, -vv, -vvv, -v4, -v5, etc.
                  Various levels of logging (-nocolor for monochrome)
  -V, --version   Print GNU-style version info
  -h, --help      Print help message
  -seed &lt;n&gt;       Seed random number generator (mt19937; default seed 5489)

REFERENCES

The reconstruction method uses phylogenetic transducers, as described in:
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

</code></pre>
