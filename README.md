# Indel Historian
Reconstruction of phylogenetic insertion-deletion histories using the transducer method
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
  -auto &lt;file&gt;    Auto-detect file format and purpose
  -model &lt;file&gt;   Specify substitution & indel model (JSON)
  -seqs &lt;file&gt;    Specify unaligned sequences (FASTA)
  -guide &lt;file&gt;   Specify guide alignment (gapped FASTA)
  -tree &lt;file&gt;    Specify phylogeny (New Hampshire)
  -nexus &lt;file&gt;, -stockholm &lt;file&gt;
                  Specify phylogeny & guide alignment together

  -saveguide &lt;file&gt;, -savemodel &lt;file&gt;
                  Save guide alignment/model to files

  -output (nexus|fasta|stockholm)
                  Specify output format (default is Stockholm)

Reconstruction algorithm options:
  -ancseq         Predict ancestral sequences (default is to leave them as *'s)
  -band &lt;n&gt;       Size of band around guide alignment (default 10)
  -noband         Turn off band, ignore guide alignment
  -minpost &lt;p&gt;    Posterior prob. threshold for profile states (default .1)
  -states &lt;n&gt;     Limit max number of states per profile

Guide alignment construction options:
  -allvsall       Try all pairwise alignments, not just a random spanning graph
  -kmatch &lt;k&gt;     Length of kmers for pre-filtering heuristic (default 6)
  -kmatchn &lt;n&gt;    Threshold# of kmer matches to seed a diagonal
  -kmatchband &lt;n&gt; Size of DP band around kmer-matching diagonals (default 64)
  -kmatchmb &lt;M&gt;   Set kmer threshold to use M megabytes of memory
  -kmatchmax      Set kmer threshold to use all available memory (default)
  -kmatchoff      No kmer threshold, do full DP

Model-fitting and event-counting options:
  -recon &lt;file&gt;, -nexusrecon&lt;file&gt;
                  Use precomputed reconstruction (FASTA/NEXUS, respectively)
  -mininc &lt;n&gt;     EM convergence threshold as relative log-likelihood increase
                    (default is .01)
  -maxiter &lt;n&gt;    Max number of EM iterations (default 100)
  -nolaplace      Do not add Laplace +1 pseudocounts during model-fitting

General options:
  -verbose, -vv, -vvv, -v4, -v5, etc.
                   Various levels of logging (-nocolor for monochrome)
  -V, --version   Print GNU-style version info
  -h, --help      Print help message
  -seed &lt;n&gt;       Seed random number generator

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
