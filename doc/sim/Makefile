
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
gp120guide.isg:
	bin/add-indelign-params.js 500 30 .3 $(GP120TREE) >$@

# General recipe for indel-seq-gen usage (and PRANK usage) is from Text S1 of Westesson et al, 2012
# http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0034572#s1
GUIDETREE := gp120guide.isg

# SUBRATESCALE is 1/1.2 to compensate for indel-seq-gen's natural substitution rate of 1.2
SUBRATESCALE := .83333

gp120sim/%.seq:
	cat $(GUIDETREE) | bin/indel-seq-gen --rng_seed `perl -e 'srand($*);print int(rand()*0xffffffff)'` -m JTT -u xia --num_gamma_cats 3 -a 0.372 --branch_scale $(SUBRATESCALE) --outfile gp120sim/$* --quiet --outfile_format f -s 10000 --write_anc

gp120sim/%.ma.json: gp120sim/%.seq
	bin/historian fit -fixsubrates -band 1 gp120sim/$*.ma >$@

gp120sim/%.hist.json: gp120sim/%.seq
	bin/historian fit -fixsubrates gp120sim/$*.seq >$@

gp120sim/%.muscle: gp120sim/%.seq
	bin/muscle -in gp120sim/$*.seq >$@
	bin/historian fit -fixsubrates -band 1 $@ >$@.json

gp120sim/%.probcons: gp120sim/%.seq
	bin/probcons gp120sim/$*.seq >$@
	bin/historian fit -fixsubrates -band 1 $@ >$@.json

gp120sim/%.prank: gp120sim/%.seq
	bin/prank -d=gp120sim/$*.seq -noxml -realbranches -writeanc -o=gp120sim/$*.prankout +F
	cp gp120sim/$*.prankout.1.fas $@
	bin/historian fit -fixsubrates -band 1 $@ >$@.json

gp120sim/%.summary: gp120sim/%.ma.json gp120sim/%.hist.json gp120sim/%.muscle gp120sim/%.probcons gp120sim/%.prank
	node -e 'fs=require("fs");process.argv.slice(1).forEach(function(arg){json=JSON.parse(fs.readFileSync(arg));console.log([arg,json.insrate,json.delrate].join(" "))})' gp120sim/$*.ma.json gp120sim/$*.hist.json gp120sim/$*.muscle.json gp120sim/$*.probcons.json gp120sim/$*.prank.json >$@


# RMS error = sqrt(sum (r-1)^2) = sum (r^2 - 2r + 1)
# where r = estimated_rate / true_rate
gp120sim/1-%.summary:
	node -e 'fs=require("fs");p=process.argv[1];i=d=-Math.log(1-p);function mean_rmserr(m1,m2){m=m1/$*;return [m,Math.sqrt(m2/$*-2*m+1)]}process.argv.slice(3).forEach(function(method){m1i=m2i=m1d=m2d=0;for(n=1;n<=$*;++n){json=JSON.parse(fs.readFileSync(process.argv[2]+n+"."+method+".json"));ir=json.insrate/i;dr=json.delrate/d;m1i+=ir;m2i+=ir*ir;m1d+=dr;m2d+=dr*dr}console.log([method].concat(mean_rmserr(m1i,m2i)).concat(mean_rmserr(m1d,m2d)).join(" "))})' 0.3 gp120sim/ ma hist muscle probcons prank >$@
	cp $@ gp120.1-100.summary

gp120sim-%:
	make `perl -e 'print map "gp120sim/$$_.summary ", 1..$*'`


# Simulation #2: DNA
gp120dna: gp120/dna/3.summary gp120/dna/4.summary gp120/dna/5.summary gp120/dna/6.summary gp120/dna/7.summary

gp120/dna/%.summary:
	test -e gp120/dna/$* || mkdir -p gp120/dna/$*
	bin/dna-sim.js $* gp120/dna/$* ../../../Makefile .3 100 $(GP120TREE)

# the following targets will be made from within the subdirectory gp120/dna/$logLength
dna%.seq:
	cat params.isg | ../../../bin/indel-seq-gen --rng_seed `perl -e 'srand($*);print int(rand()*0xffffffff)'` -m GTR -u ran --branch_scale $(SUBRATESCALE) --outfile dna$* --quiet --outfile_format f -s 10000 --write_anc

dna%.hist.json: dna%.seq
	../../../bin/historian fit -fixsubrates dna$*.seq -v2 >$@

# keep intermediates
.SECONDARY:
