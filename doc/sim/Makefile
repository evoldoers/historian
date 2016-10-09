
MAKEDIR: $(dir $(abspath $(lastword $(MAKEFILE_LIST))))

progs: bin/rose bin/prank bin/muscle bin/probcons bin/historian bin/qscore

bin/rose:
	test -e bin || mkdir bin
	curl -O https://bibiserv2.cebitec.uni-bielefeld.de/applications/rose/resources/downloads/rose-1.3.1-src.tar.gz
	tar xvzf rose-1.3.1-src.tar.gz
	cd rose-1.3.1-src; ./configure --prefix=$(MAKEDIR); make; make install

bin/prank:
	test -e bin || mkdir bin
	curl -O http://www.ebi.ac.uk/goldman-srv/prank/src/prank/prank.osx_1058.100701.tgz
	tar xvzf prank.osx_1058.100701.tgz

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
