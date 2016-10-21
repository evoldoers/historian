
# Installation

MAKEDIR := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))

progs: bin/indel-seq-gen bin/rose bin/prank bin/muscle bin/probcons bin/historian bin/qscore
	rm *.tar.gz *.tgz

bin/indel-seq-gen:
	curl -O http://bioinfolab.unl.edu/~cstrope/iSG/indel-seq-gen-2.1.03.tar.gz
	tar xvzf indel-seq-gen-2.1.03.tar.gz
	perl -pi.bak -e 'print"#include <unistd.h>\n"unless$$n++' indel-seq-gen-2.1.03/src/main.cpp
	perl -pi.bak -e 's/next/nxt/' indel-seq-gen-2.1.03/src/random.cpp
	cd indel-seq-gen-2.1.03; ./configure --prefix=$(MAKEDIR); make; make install

bin/rose:
	curl -O https://bibiserv2.cebitec.uni-bielefeld.de/applications/rose/resources/downloads/rose-1.3.1-src.tar.gz
	tar xvzf rose-1.3.1-src.tar.gz
	cd rose-1.3.1; ./configure --prefix=$(MAKEDIR); make; make install

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
# These are then used to make gp120guide.isg
#  500 is the length of the root domain (GP120)
#   30 is the max indel length (by inspection of gp120.recon.fa)
# All these parameters are rounded to 1 significant figure
gp120.json:
	bin/historian fit ../../data/gp120.fa >$@

gp120.recon.fa: gp120.json
	bin/historian ../../data/gp120.fa -model gp120.json >$@

GP120TREE := ../../data/gp120.tree.nh
SYM8TREE  := sym8tree.nh

GUIDETREE := tree.nh

# General recipe for indel-seq-gen usage (and PRANK usage) is from Text S1 of Westesson et al, 2012
# http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0034572#s1
# but simplified a little (no substitution rate scaling)
ISGPARAMS := guide.isg

# INDELRATE and INDELPROB
# default rate is by inspection of gp120.json
GP120JSON = gp120.json
INDELRATE ?= $(shell node -e 'fs=require("fs");json=JSON.parse(fs.readFileSync("$(GP120JSON)"));process.stdout.write(((json.insrate + json.delrate) / 2).toPrecision(2))')
INDELPROB := $(shell perl -e 'print 1-exp(-$(INDELRATE))')

RATEDIR := gp120sim/rate$(INDELRATE)

indelprob-test:
	echo $(INDELRATE)
	echo $(INDELPROB)

# sequence length & max indel length
SEQLEN := 500
MAXINDELLEN := 30

$(RATEDIR)/tree%/$(GUIDETREE):
	test -e $(dir $@) || mkdir -p $(dir $@)
	bin/scale-tree.js $(GP120TREE) $* >$(dir $@)/tree.nh

$(RATEDIR)/sym%/$(GUIDETREE):
	test -e $(dir $@) || mkdir -p $(dir $@)
	bin/scale-tree.js $(SYM8TREE) $* >$(dir $@)/tree.nh

$(RATEDIR)/%/$(ISGPARAMS): $(RATEDIR)/%/$(GUIDETREE)
	bin/add-indelign-params.js $(SEQLEN) $(MAXINDELLEN) $(INDELPROB) $(dir $@)/tree.nh $* >$@

# Simulation
%.seq:
	make $(dir $@)$(ISGPARAMS)
	cat $(dir $@)$(ISGPARAMS) | bin/indel-seq-gen --rng_seed `perl -e 'srand($(notdir $*));print int(rand()*0xffffffff)'` -m JTT -u xia --outfile $* --quiet --outfile_format f -s 10000 --write_anc

%.leaf: %.seq
	cat $< | perl -e 'while(<>){if(/^>.*/){$$n=$$_}print if$$n && !($$n=~/^>\d+\s*$$/)}' >$@

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

%.histref.stk: %.leaf
	bin/historian recon -refine $*.leaf -v3 >$@

%.histref.json: %.histref.stk
	bin/historian fit -fixsubrates -stockrecon $*.histref.stk -v3 >$@

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

%.summary: %.ma.json %.ma.tree.json %.hist.json %.histrec.json %.histref.json %.histslow.json %.muscle %.probcons %.prank
	node -e 'fs=require("fs");process.argv.slice(1).forEach(function(arg){json=JSON.parse(fs.readFileSync(arg));console.log([arg,json.insrate,json.delrate].join(" "))})' $*.ma.json $*.ma.tree.json $*.hist.json $*.histrec.json $*.histref.json $*.histslow.json $*.muscle.json $*.probcons.json $*.prank.json >$@

# RMS error = sqrt(sum (r-1)^2) = sum (r^2 - 2r + 1)
# where r = estimated_rate / true_rate
REPS ?= 100
%/1-$(REPS).summary:
	node -e 'fs=require("fs");p=process.argv[1];i=d=-Math.log(1-p);function mean_rmserr(m1,m2){m=m1/$(REPS);return [m,Math.sqrt(m2/$(REPS)-2*m+1)]}process.argv.slice(3).forEach(function(method){m1i=m2i=m1d=m2d=0;for(n=1;n<=$(REPS);++n){json=JSON.parse(fs.readFileSync(process.argv[2]+n+"."+method+".json"));ir=json.insrate/i;dr=json.delrate/d;m1i+=ir;m2i+=ir*ir;m1d+=dr;m2d+=dr*dr}console.log([method].concat(mean_rmserr(m1i,m2i)).concat(mean_rmserr(m1d,m2d)).join(" "))})' $(INDELPROB) $(dir $@) ma ma.tree hist histrec histref histslow muscle probcons prank >$@

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
