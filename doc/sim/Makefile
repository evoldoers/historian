
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
#   .3 is the indel opening probability = 1 - exp(-(insRate+delRate)/2)
# All these parameters are rounded to 1 significant figure
gp120.json:
	bin/historian fit ../../data/gp120.fa >$@

gp120.recon.fa: gp120.json
	bin/historian ../../data/gp120.fa -model gp120.json >$@

GP120TREE := ../../data/gp120.tree.nh
SYM8TREE  := sym8tree.nh

# General recipe for indel-seq-gen usage (and PRANK usage) is from Text S1 of Westesson et al, 2012
# http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0034572#s1
GUIDETREE := guide.isg

# SUBRATESCALE is 1/1.2 to compensate for indel-seq-gen's natural substitution rate of 1.2
SUBRATESCALE := .83333

gp120sim/tree%/$(GUIDETREE):
	test -e gp120sim/tree$* || mkdir -p gp120sim/tree$*
	bin/add-indelign-params.js 500 30 .3 $(GP120TREE) $* >$@

gp120sim/sym%/$(GUIDETREE):
	test -e gp120sim/sym$* || mkdir -p gp120sim/sym$*
	bin/add-indelign-params.js 500 30 .3 $(SYM8TREE) $* >$@

%.seq:
	make $(dir $@)$(GUIDETREE)
	cat $(dir $@)$(GUIDETREE) | bin/indel-seq-gen --rng_seed `perl -e 'srand($(notdir $*));print int(rand()*0xffffffff)'` -m JTT -u xia --num_gamma_cats 3 -a 0.372 --branch_scale $(SUBRATESCALE) --outfile $* --quiet --outfile_format f -s 10000 --write_anc

%.ma.json: %.seq
	bin/historian fit -fixsubrates -band 1 $*.ma >$@

%.hist.json: %.seq
	bin/historian fit -fixsubrates $*.seq -v3 >$@

%.histrec.stk: %.seq
	bin/historian recon $*.seq -v3 >$@

%.histrec.json: %.histrec.stk
	bin/historian fit -fixsubrates -stockrecon $*.histrec.stk -v3 >$@

%.muscle: %.seq
	bin/muscle -in $*.seq >$@
	bin/historian fit -fixsubrates -band 1 $@ >$@.json

%.probcons: %.seq
	bin/probcons $*.seq >$@
	bin/historian fit -fixsubrates -band 1 $@ >$@.json

%.prank: %.seq
	bin/prank -d=$*.seq -noxml -realbranches -writeanc -o=$*.prankout +F
	cp $*.prankout.1.fas $@
	bin/historian fit -fixsubrates -band 1 $@ >$@.json

%.summary: %.ma.json %.hist.json %.histrec.json %.muscle %.probcons %.prank
	node -e 'fs=require("fs");process.argv.slice(1).forEach(function(arg){json=JSON.parse(fs.readFileSync(arg));console.log([arg,json.insrate,json.delrate].join(" "))})' $*.ma.json $*.hist.json $*.muscle.json $*.probcons.json $*.prank.json >$@


# RMS error = sqrt(sum (r-1)^2) = sum (r^2 - 2r + 1)
# where r = estimated_rate / true_rate
REPS := 100
%/1-$(REPS).summary:
	node -e 'fs=require("fs");p=process.argv[1];i=d=-Math.log(1-p);function mean_rmserr(m1,m2){m=m1/$(REPS);return [m,Math.sqrt(m2/$(REPS)-2*m+1)]}process.argv.slice(3).forEach(function(method){m1i=m2i=m1d=m2d=0;for(n=1;n<=$(REPS);++n){json=JSON.parse(fs.readFileSync(process.argv[2]+n+"."+method+".json"));ir=json.insrate/i;dr=json.delrate/d;m1i+=ir;m2i+=ir*ir;m1d+=dr;m2d+=dr*dr}console.log([method].concat(mean_rmserr(m1i,m2i)).concat(mean_rmserr(m1d,m2d)).join(" "))})' 0.3 $(dir $@) ma hist histrec muscle probcons prank >$@

gp120sim-tree%:
	make `perl -e 'print map "gp120sim/tree$*/$$_.summary ", 1..$(REPS)'`

gp120sim-sym%:
	make `perl -e 'print map "gp120sim/sym$*/$$_.summary ", 1..$(REPS)'`


# Simulation #2: DNA
gp120dna: gp120sim/dna/3/dna.hist.stk gp120sim/dna/4/dna.hist.stk gp120sim/dna/5/dna.hist.stk gp120sim/dna/6/dna.hist.stk gp120sim/dna/7/dna.hist.stk

gp120sim/dna/%/$(GUIDETREE):
	test -e $(dir $@) || mkdir -p $(dir $@)
	bin/add-indelign-params.js `perl -e 'print 10 ** $*'`  `perl -e 'print 2 ** ($* + 2)'` .3 $(SYM8TREE) >$@
	cat $@

%/dna: %/$(GUIDETREE)
	cat $(dir $@)$(GUIDETREE) | bin/indel-seq-gen --rng_seed `perl -e 'srand($(notdir $*));print int(rand()*0xffffffff)'` -m GTR -u ran --branch_scale $(SUBRATESCALE) --outfile $*/dna --quiet --outfile_format f -s 10000 --write_anc
	mv $*/dna.seq $@

%/dna.hist.stk: %/dna
	bin/historian recon $*/dna -fast -v2 >$@

# keep intermediates
.SECONDARY:
