
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
	test -e bin || mkdir bin
	curl -O https://bibiserv2.cebitec.uni-bielefeld.de/applications/rose/resources/downloads/rose-1.3.1-src.tar.gz
	tar xvzf rose-1.3.1-src.tar.gz
	cd rose-1.3.1; ./configure --prefix=$(MAKEDIR); make; make install

bin/prank:
	test -e bin || mkdir bin
	curl -O http://www.ebi.ac.uk/goldman-srv/prank/src/prank/prank.osx_1058.100701.tgz
	tar xvzf prank.osx_1058.100701.tgz
	mv prank $@

bin/muscle:
	test -e bin || mkdir bin
	curl -O http://www.drive5.com/muscle/downloads3.8.31/muscle3.8.31_i86darwin64.tar.gz
	tar xvzf muscle3.8.31_i86darwin64.tar.gz
	mv muscle3.8.31_i86darwin64 $@

bin/probcons:
	test -e bin || mkdir bin
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


# Simulation

# General formula is from Text S1 of Westesson et al, 2012
# http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0034572#s1
GUIDETREE := gp120guide.isg

# SUBRATESCALE is 1/1.2 to compensate for indel-seq-gen's natural substitution rate of 1.2
SUBRATESCALE := .83333

gp120sim:
	cat $(GUIDETREE) | bin/indel-seq-gen -m JTT -u xia --num_gamma_cats 3 -a 0.372 --branch_scale $(SUBRATESCALE) --outfile $@ --quiet --outfile_format f -s 10000 --write_anc
