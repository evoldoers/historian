# Indel Historian
Reconstruction of phylogenetic insertion-deletion histories

<pre><code>
Usage: historian {align,help,version} [options]

Reconstruction:

  historian align seqs.fa [-tree tree.nh] [-model model.json] &gt;alignment.fa

Reconstruction options:
   -seqs &lt;file&gt;    Specify unaligned sequences (FASTA)
   -guide &lt;file&gt;   Specify guide alignment (gapped FASTA)
   -tree &lt;file&gt;    Specify phylogeny (New Hampshire)
   -model &lt;file&gt;   Specify substitution & indel model (JSON)

   -saveseqs &lt;file&gt;, -saveguide &lt;file&gt;, -savetree &lt;file&gt;, -savemodel &lt;file&gt;
                   Save intermediate analysis datasets to specified files

   -band &lt;n&gt;       Size of band around guide alignment
   -noband         Turn off band, ignore guide alignment
   -samples &lt;n&gt;    Max number of alignments to sample when building profiles
   -states &lt;n&gt;     Max number of states allowed in any single profile
   -seed &lt;n&gt;       Seed random number generator

Guide alignment options:
   -kmatch &lt;k&gt;     Length of kmers for pre-filtering heuristic (default 6)
   -kmatchn &lt;n&gt;    Threshold# of kmer matches to seed a diagonal
   -kmatchband &lt;n&gt; Size of DP band around kmer-matching diagonals (default 64)
   -kmatchmb &lt;M&gt;   Set kmer threshold to use M megabytes of memory
   -kmatchmax      Set kmer threshold to use all available memory (default)
   -kmatchoff      No kmer threshold, do full DP

</code></pre>
