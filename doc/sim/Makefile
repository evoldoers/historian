
# Installation

MAKEDIR := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
HISTDIR := $(MAKEDIR)/../..

progs: bin/prank bin/muscle bin/probcons bin/historian bin/qscore
	rm *.tar.gz *.tgz

bin/prank:
	curl -O http://www.ebi.ac.uk/goldman-srv/prank/src/prank/prank.osx_1058.100701.tgz
	tar xvzf prank.osx_1058.100701.tgz
	mv prank $@

bin/muscle:
	curl -O http://www.drive5.com/muscle/downloads3.8.31/muscle3.8.31_i86darwin64.tar.gz
	tar xvzf muscle3.8.31_i86darwin64.tar.gz
	mv muscle3.8.31_i86darwin64 $@

bin/probcons:
	curl -O http://probcons.stanford.edu/probcons_v1_12.tar.gz
	tar xvzf probcons_v1_12.tar.gz
	cd probcons; make
	mv probcons/probcons $@

bin/historian:
	cd $(HISTDIR); make historian
	cp $(HISTDIR)/bin/historian $@

bin/qscore:
	test -e qscore || mkdir qscore
	cd qscore; curl -O http://www.drive5.com/qscore/qscore_src.tar.gz
	cd qscore; tar xvzf qscore_src.tar.gz
	cd qscore; make
	mv qscore/qscore $@


# Simulation #1: GP120 protein-like

# First estimate rough parameters of indel evolution on GP120 tree
gp120.json:
	bin/historian -fast fit $(HISTDIR)/data/gp120.fa >$@

gp120.recon.fa: gp120.json
	bin/historian -fast $(HISTDIR)/data/gp120.fa -model gp120.json >$@

# Simulation model parameters
PRESET ?= dayhoff
GP120JSON = gp120.json
INDELRATE ?= $(shell node -e 'fs=require("fs");json=JSON.parse(fs.readFileSync("$(GP120JSON)"));process.stdout.write(((json.insrate + json.delrate) / 2).toPrecision(2))')
INDELEXTPROB ?= $(shell node -e 'fs=require("fs");json=JSON.parse(fs.readFileSync("$(GP120JSON)"));i=json.insextprob;d=json.delextprob;process.stdout.write((2/(1/i+1/d)).toPrecision(2))')
SEQLEN := 500

SIMPARAMS := simparams.txt

RATEDIR := gp120sim/rate$(INDELRATE)
GP120TREE := $(HISTDIR)/data/gp120.tree.nh
SYM8TREE  := sym8tree.nh
GUIDETREE := tree.nh

$(RATEDIR)/tree%/$(GUIDETREE):
	test -e $(dir $@) || mkdir -p $(dir $@)
	bin/scale-tree.js $(GP120TREE) $* >$(dir $@)/tree.nh

$(RATEDIR)/sym%/$(GUIDETREE):
	test -e $(dir $@) || mkdir -p $(dir $@)
	bin/scale-tree.js $(SYM8TREE) $* >$(dir $@)/tree.nh

$(RATEDIR)/%/$(SIMPARAMS): $(RATEDIR)/%/$(GUIDETREE)
	echo "-preset $(PRESET) -gaprate $(INDELRATE) -gapextprob $(INDELEXTPROB)" >$@

# Simulation
$(RATEDIR)/%.ma.stk:
	make $(dir $@)$(SIMPARAMS) $(dir $@)$(GUIDETREE)
	bin/historian generate -seed `perl -e 'srand($(notdir $*));print int(rand()*0xffffffff)'` `cat $(dir $@)$(SIMPARAMS)` $(dir $@)$(GUIDETREE) -rootlen $(SEQLEN) -savemodel $(RATEDIR)/$*.sim.json >$@

%.ma: %.ma.stk
	$(HISTDIR)/perl/stockholm2fasta.pl -g $< >$@

%.seq: %.ma.stk
	$(HISTDIR)/perl/stockholm2fasta.pl $< >$@

%.leaf: %.seq
	cat $< | perl -e 'while(<>){if(/^>.*/){$$n=$$_}print if$$n && !($$n=~/^>node\d+\s*$$/)}' >$@

%.ma.json: %.leaf
	bin/historian fit -fast -fixsubrates -band 1 $*.ma >$@

%.ma.tree.json: %.leaf
	bin/historian fit -fast -fixsubrates -band 1 $*.ma -tree $(dir $@)/tree.nh >$@

%.hist.json: %.leaf
	bin/historian fit -fast -fixsubrates $*.leaf -v3 >$@

%.histslow.json: %.leaf
	bin/historian fit -fixsubrates -profminpost .001 -profmaxstates 12000 $*.leaf -v3 >$@

%.histrec.stk: %.leaf
	bin/historian recon -fast -norefine $*.leaf -v3 >$@

%.histrec.json: %.histrec.stk
	bin/historian fit -fast -fixsubrates -stockrecon $*.histrec.stk -v3 >$@

%.histref.stk: %.leaf
	bin/historian recon -fast -refine $*.leaf -v3 >$@

%.histref.json: %.histref.stk
	bin/historian fit -fast -fixsubrates -stockrecon $*.histref.stk -v3 >$@

%.histfit.stk: %.hist.json
	bin/historian recon -fast -model $*.hist.json $*.leaf -v3 >$@

%.histfit.json: %.histfit.stk
	bin/historian fit -fast -fixsubrates -model $*.hist.json -stockrecon $*.histfit.stk -v3 >$@

%.histfit.fa: %.histfit.stk
	$(HISTDIR)/perl/stockholm2fasta.pl -g $< >$@

%.histprk.stk: %.leaf
	bin/historian recon -fast -preset wag -gaprate .005 -gapextprob .5 -keepgapsopen $*.leaf -v3 >$@

%.histprk.json: %.histprk.stk
	bin/historian fit -fast -preset wag -gaprate .005 -gapextprob .5 -fixsubrates -stockrecon $< -v3 >$@

%.histrand.stk: %.hist.json
	bin/historian recon -fast -model $*.hist.json -profsamples 100 $*.leaf -v3 >$@

%.histrand.json: %.histrand.stk
	bin/historian fit -fast -fixsubrates -model $*.hist.json -stockrecon $< -v3 >$@

%.muscle.json: %.leaf
	bin/muscle -in $*.leaf >$*.muscle
	bin/historian fit -fast -fixsubrates -band 1 $*.muscle >$@

%.probcons.json: %.leaf
	bin/probcons $*.leaf >$*.probcons
	bin/historian fit -fast -fixsubrates -band 1 $*.probcons >$@

%.prank.json: %.leaf
	bin/prank -d=$*.leaf -noxml -realbranches -writeanc -o=$*.prankout +F
	cp $*.prankout.1.fas $*.prank
	bin/historian fit -fast -fixsubrates -band 1 $*.prank >$@

ALLMETHODS := ma ma.tree hist histrec histref histslow histfit histprk histrand muscle probcons prank
ALLJSON := $(addprefix %.,$(addsuffix .json,$(ALLMETHODS)))

%.summary: $(ALLJSON)
	node -e 'fs=require("fs");process.argv.slice(1).forEach(function(arg){json=JSON.parse(fs.readFileSync(arg));console.log([arg,json.insrate,json.delrate].join(" "))})' $*.ma.json $*.ma.tree.json $*.hist.json $*.histrec.json $*.histref.json $*.histslow.json $*.histfit.json $*.muscle.json $*.probcons.json $*.prank.json >$@

# RMS error = sqrt((1/N) sum (r-1)^2) = sqrt(sum ((r^2)/N - 2r/N + 1))
# where r = estimated_rate / true_rate
REPS ?= 100
%/1-$(REPS).summary:
	bin/summary.js "$(ALLMETHODS)" $* >$@

gp120sim/all.1-$(REPS).summary:
	bin/summary.js "$(ALLMETHODS)" $(wildcard gp120sim/rate*/*) >$@

gp120sim/alltree.1-$(REPS).summary:
	bin/summary.js "$(ALLMETHODS)" $(wildcard gp120sim/rate*/tree*) >$@

gp120sim/rate%/alltree.1-$(REPS).summary:
	bin/summary.js "$(ALLMETHODS)" $(wildcard gp120sim/rate$*/tree*) >$@

gp120sim/rate%/allsym.1-$(REPS).summary:
	bin/summary.js "$(ALLMETHODS)" $(wildcard gp120sim/rate$*/sym*) >$@

gp120sim/allsym.1-$(REPS).summary:
	bin/summary.js "$(ALLMETHODS)" $* $(wildcard gp120sim/rate*/sym*) >$@

gp120sim/results.1-$(REPS).dat:
	bin/results.js .1 2 false "$(ALLMETHODS)" $(wildcard gp120sim/rate*/*) >$@

gp120sim/together.1-$(REPS).dat:
	bin/results.js .1 2 true "$(ALLMETHODS)" $(wildcard gp120sim/rate*/*) >$@

allsummary: $(addsuffix /1-$(REPS).summary,$(wildcard gp120sim/rate*/*))

gp120sim-tree%:
	make INDELRATE=$(INDELRATE) `perl -e 'print map "gp120sim/rate$(INDELRATE)/tree$*/$$_.summary ", 1..$(REPS)'`

gp120sim-sym%:
	make INDELRATE=$(INDELRATE) `perl -e 'print map "gp120sim/rate$(INDELRATE)/sym$*/$$_.summary ", 1..$(REPS)'`

# Top-level targets for simulation #1
sim1:
	make INDELRATE=0.005 gp120sim-tree5
	make INDELRATE=0.005 gp120sim-sym5
	make INDELRATE=0.01 gp120sim-tree5
	make INDELRATE=0.01 gp120sim-sym5
	make INDELRATE=0.028 gp120sim-tree1
	make INDELRATE=0.028 gp120sim-sym1
	make INDELRATE=0.028 gp120sim-tree2
	make INDELRATE=0.028 gp120sim-sym2
	make INDELRATE=0.028 gp120sim-tree5
	make INDELRATE=0.028 gp120sim-sym5
	make INDELRATE=0.056 gp120sim-tree2
	make INDELRATE=0.056 gp120sim-sym2
	make INDELRATE=0.056 gp120sim-tree5
	make INDELRATE=0.056 gp120sim-sym5
	make results.pdf

results.pdf:
	R -f bin/plotresults.R

# keep intermediates
.SECONDARY:
