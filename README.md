# Indel Historian
Reconstruction of phylogenetic insertion-deletion histories using the transducer method
(see [Westesson et al, 2012](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0034572) for an evaluation and brief description of the method, or [this arXiv report](http://arxiv.org/abs/1103.4347) for a tutorial introduction).

<pre><code>
Usage: historian {reconstruct,count,post,sum,fit,help,version} [options]

COMMANDS

Reconstruction:
  historian reconstruct seqs.fa [-tree tree.nh] &gt;reconstruction.fa
  historian reconstruct -guide guide.fa [-tree tree.nh] &gt;reconstruction.fa

Event counting:
  historian count reconstruction.fa tree.nh [-model model.json] &gt;counts.json
  historian post seqs.fa [-tree tree.nh] [-model model.json] &gt;counts.json

The former (count) obtains a point estimate from a single reconstruction.
The latter (post) averages over (a subset of) the posterior distribution.

Model fitting:
  historian sum counts.json pseudocounts.json [morecounts.json...] &gt;sum.json
  historian fit counts.json &gt;newmodel.json

All commands can be abbreviated to single letters, like so:
  historian r seqs.fa &gt;reconstruction.fa
  historian p seqs.fa &gt;counts.json
  historian f counts.json &gt;model.json
(etc.)

OPTIONS

Reconstruction file I/O options:
  -seqs &lt;file&gt;    Specify unaligned sequences (FASTA)
  -guide &lt;file&gt;   Specify guide alignment (gapped FASTA)
  -tree &lt;file&gt;    Specify phylogeny (New Hampshire)
  -model &lt;file&gt;   Specify substitution & indel model (JSON)

  -saveseqs &lt;file&gt;, -saveguide &lt;file&gt;, -savetree &lt;file&gt;, -savemodel &lt;file&gt;
                  Save various intermediate analysis results to files

Reconstruction algorithm options:
  -ancseq         Predict ancestral sequences (default is to leave them as *'s)
  -band &lt;n&gt;       Size of band around guide alignment (default 10)
  -noband         Turn off band, ignore guide alignment
  -minpost &lt;p&gt;    Posterior prob. threshold for profile states (default .1)
  -states &lt;n&gt;     Limit max number of states per profile

Guide alignment construction options:
  -kmatch &lt;k&gt;     Length of kmers for pre-filtering heuristic (default 6)
  -kmatchn &lt;n&gt;    Threshold# of kmer matches to seed a diagonal
  -kmatchband &lt;n&gt; Size of DP band around kmer-matching diagonals (default 64)
  -kmatchmb &lt;M&gt;   Set kmer threshold to use M megabytes of memory
  -kmatchmax      Set kmer threshold to use all available memory (default)
  -kmatchoff      No kmer threshold, do full DP

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
