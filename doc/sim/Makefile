
# Installation

MAKEDIR := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))

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
	cd ../..; make historian
	cp ../../bin/historian $@

bin/qscore:
	test -e qscore || mkdir qscore
	cd qscore; curl -O http://www.drive5.com/qscore/qscore_src.tar.gz
	cd qscore; tar xvzf qscore_src.tar.gz
	cd qscore; make
	mv qscore/qscore $@


# Simulation #1: GP120 protein-like

# First estimate rough parameters of indel evolution on GP120 tree
gp120.json:
	bin/historian fit ../../data/gp120.fa >$@

gp120.recon.fa: gp120.json
	bin/historian ../../data/gp120.fa -model gp120.json >$@

# Simulation model parameters
PRESET ?= dayhoff
GP120JSON = gp120.json
INDELRATE ?= $(shell node -e 'fs=require("fs");json=JSON.parse(fs.readFileSync("$(GP120JSON)"));process.stdout.write(((json.insrate + json.delrate) / 2).toPrecision(2))')
INDELEXTPROB ?= $(shell node -e 'fs=require("fs");json=JSON.parse(fs.readFileSync("$(GP120JSON)"));i=json.insextprob;d=json.delextprob;process.stdout.write((2/(1/i+1/d)).toPrecision(2))')
SEQLEN := 500

SIMPARAMS := simparams.txt

RATEDIR := gp120sim/rate$(INDELRATE)
GP120TREE := ../../data/gp120.tree.nh
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
	../../perl/stockholm2fasta.pl -g $< >$@

%.seq: %.ma.stk
	../../perl/stockholm2fasta.pl $< >$@

%.leaf: %.seq
	cat $< | perl -e 'while(<>){if(/^>.*/){$$n=$$_}print if$$n && !($$n=~/^>node\d+\s*$$/)}' >$@

%.ma.json: %.leaf
	bin/historian fit -fixsubrates -band 1 $*.ma >$@

%.ma.tree.json: %.leaf
	bin/historian fit -fixsubrates -band 1 $*.ma -tree $(dir $@)/tree.nh >$@

%.hist.json: %.leaf
	bin/historian fit -fixsubrates $*.leaf -v3 >$@

%.histslow.json: %.leaf
	bin/historian fit -fixsubrates -profminpost .001 -profmaxstates 12000 $*.leaf -v3 >$@

%.histrec.stk: %.leaf
	bin/historian recon -norefine $*.leaf -v3 >$@

%.histrec.json: %.histrec.stk
	bin/historian fit -fixsubrates -stockrecon $*.histrec.stk -v3 >$@

%.histref.stk: %.leaf
	bin/historian recon -refine $*.leaf -v3 >$@

%.histref.json: %.histref.stk
	bin/historian fit -fixsubrates -stockrecon $*.histref.stk -v3 >$@

%.histfit.stk: %.hist.json
	bin/historian recon -model $*.hist.json $*.leaf -v3 >$@

%.histfit.json: %.histfit.stk
	bin/historian fit -fixsubrates -model $*.hist.json -stockrecon $*.histfit.stk -v3 >$@

%.histfit.fa: %.histfit.stk
	../../perl/stockholm2fasta.pl -g $< >$@

%.hist2.json: %.histfit.fa
	bin/historian fit -fixsubrates -model $*.hist.json -guide $*.histfit.fa -v3 >$@

%.histprk.stk: %.leaf
	bin/historian recon -preset wag -gaprate .005 -gapextprob .5 -keepgapsopen $*.leaf -v3 >$@

%.histprk.json: %.histprk.stk
	bin/historian fit -preset wag -gaprate .005 -gapextprob .5 -fixsubrates -stockrecon $< -v3 >$@

%.histmcmc.stk: %.histfit.stk %.hist.json
	bin/historian mcmc -fixtree -model $*.hist.json $*.histfit.stk -v4 >$@

%.histmcmc.json: %.histmcmc.stk
	bin/historian fit -stockrecon $< -v3 >$@

%.histrand.stk: %.hist.json
	bin/historian recon -model $*.hist.json -profsamples 100 $*.leaf -v3 >$@

%.histrand.json: %.histrand.stk
	bin/historian fit -fixsubrates -model $*.hist.json -stockrecon $< -v3 >$@

%.muscle: %.leaf
	bin/muscle -in $*.leaf >$@
	bin/historian fit -fixsubrates -band 1 $@ >$@.json

%.probcons: %.leaf
	bin/probcons $*.leaf >$@
	bin/historian fit -fixsubrates -band 1 $@ >$@.json

%.prank: %.leaf
	bin/prank -d=$*.leaf -noxml -realbranches -writeanc -o=$*.prankout +F
	cp $*.prankout.1.fas $@
	bin/historian fit -fixsubrates -band 1 $@ >$@.json

%.summary: %.ma %.ma.json %.ma.tree.json %.hist.json %.histrec.json %.histref.json %.histslow.json %.histfit.json %.hist2 %.histprk %.histmcmc %.histrand %.muscle %.probcons %.prank
	node -e 'fs=require("fs");process.argv.slice(1).forEach(function(arg){json=JSON.parse(fs.readFileSync(arg));console.log([arg,json.insrate,json.delrate].join(" "))})' $*.ma.json $*.ma.tree.json $*.hist.json $*.histrec.json $*.histref.json $*.histslow.json $*.histfit.json $*.muscle.json $*.probcons.json $*.prank.json >$@

# RMS error = sqrt((1/N) sum (r-1)^2) = sqrt(sum ((r^2)/N - 2r/N + 1))
# where r = estimated_rate / true_rate
REPS ?= 100
gp120sim/rate%/1-$(REPS).summary:
	node -e 'fs=require("fs");i=d=parseFloat("$*");dir=process.argv[1];files=process.argv.slice(2);function mean_rmserr(m1,m2){m=m1/$(REPS);return [m,Math.sqrt(m2/$(REPS)-2*m+1)]}files.forEach(function(method){m1i=m2i=m1d=m2d=0;for(n=1;n<=$(REPS);++n){json=JSON.parse(fs.readFileSync(dir+n+"."+method+".json"));ir=json.insrate/i;dr=json.delrate/d;m1i+=ir;m2i+=ir*ir;m1d+=dr;m2d+=dr*dr}console.log([method].concat(mean_rmserr(m1i,m2i)).concat(mean_rmserr(m1d,m2d)).join(" "))})' $(dir $@) ma ma.tree hist histrec histref histslow histfit hist2 histprk histmcmc histrand muscle probcons prank >$@

%/histfit:
	make `perl -e 'print map "$*/$$_.histfit.json ", 1..$(REPS)'`

%/hist2:
	make `perl -e 'print map "$*/$$_.hist2.json ", 1..$(REPS)'`

%/prk:
	make `perl -e 'print map "$*/$$_.histprk.json ", 1..$(REPS)'`

%/mcmc:
	make `perl -e 'print map "$*/$$_.histmcmc.json ", 1..$(REPS)'`

%/rand:
	make `perl -e 'print map "$*/$$_.histrand.json ", 1..$(REPS)'`

all-hist2: $(patsubst %/1-$(REPS).summary,%/hist2,$(wildcard gp120sim/*/*/1-$(REPS).summary))

all-prk: $(patsubst %/1-$(REPS).summary,%/histprk,$(wildcard gp120sim/*/*/1-$(REPS).summary))

all-mcmc: $(patsubst %/1-$(REPS).summary,%/histmcmc,$(wildcard gp120sim/*/*/1-$(REPS).summary))

all-rand: $(patsubst %/1-$(REPS).summary,%/histrand,$(wildcard gp120sim/*/*/1-$(REPS).summary))

gp120sim-tree%:
	make INDELRATE=$(INDELRATE) `perl -e 'print map "gp120sim/rate$(INDELRATE)/tree$*/$$_.summary ", 1..$(REPS)'`

gp120sim-sym%:
	make INDELRATE=$(INDELRATE) `perl -e 'print map "gp120sim/rate$(INDELRATE)/sym$*/$$_.summary ", 1..$(REPS)'`

# Top-level targets for simulation #1
sim1:
	make gp120sim-tree1
	make gp120sim-sym1
	make INDELRATE=0.014 gp120sim-tree1
	make INDELRATE=0.014 gp120sim-sym1
	make gp120sim-tree2
	make gp120sim-sym2

# keep intermediates
.SECONDARY:
