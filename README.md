# Indel Historian
Reconstruction of phylogenetic insertion-deletion histories using the transducer method
(see [here](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0034572) for an evaluation and brief description of the method, or [here](http://arxiv.org/abs/1103.4347) for a longer description).

<pre><code>
Usage: historian {align,help,version} [options]

Reconstruction:

  historian align seqs.fa [-tree tree.nh] [-model model.json] &gt;alignment.fa

File options:
   -seqs &lt;file&gt;    Specify unaligned sequences (FASTA)
   -guide &lt;file&gt;   Specify guide alignment (gapped FASTA)
   -tree &lt;file&gt;    Specify phylogeny (New Hampshire)
   -model &lt;file&gt;   Specify substitution & indel model (JSON)

   -saveseqs &lt;file&gt;, -saveguide &lt;file&gt;, -savetree &lt;file&gt;, -savemodel &lt;file&gt;
                   Save various intermediate analysis results to files

Reconstruction options:
   -band &lt;n&gt;       Size of band around guide alignment (default 10)
   -noband         Turn off band, ignore guide alignment
   -samples &lt;n&gt;    Number of sampled alignments per profile (default 100)
   -states &lt;n&gt;     Limit max number of states per profile
   -seed &lt;n&gt;       Seed random number generator

Guide alignment options:
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

The method is that of phylogenetic transducers, as described in:
 Westesson, Lunter, Paten & Holmes
 Accurate Reconstruction of Insertion-Deletion Histories by Statistical Phylogenetics
 PLoS One, DOI: 10.1371/journal.pone.0034572
 http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0034572

See also:
 Phylogenetic automata, pruning, and multiple alignment
 http://arxiv.org/abs/1103.4347
</code></pre>
