# -*- makefile-gmake -*-

.SECONDARY:

# Pseudotargets that control compilation
USING_BOOST = $(findstring boost,$(MAKECMDGOALS))
IS_DEBUG = $(findstring debug,$(MAKECMDGOALS))

# try to figure out where GSL is
# autoconf would be better but we just need a quick hack for now :)
# Thanks to Torsten Seemann for gsl-config and pkg-config formulae
GSLPREFIX = $(shell gsl-config --prefix)
ifeq (,$(wildcard $(GSLPREFIX)/include/gsl/gsl_sf.h))
GSLPREFIX = /usr
ifeq (,$(wildcard $(GSLPREFIX)/include/gsl/gsl_sf.h))
GSLPREFIX = /usr/local
endif
endif

GSLFLAGS = $(shell pkg-config --cflags gsl)
ifeq (, $(GSLFLAGS))
GSLFLAGS = -I$(GSLPREFIX)/include
endif

GSLLIBS = $(shell pkg-config --libs gsl)
ifeq (, $(GSLLIBS))
GSLLIBS = -L$(GSLPREFIX)/lib -lgsl -lgslcblas -lm
endif

# figure out whether to use Boost
# Boost is needed for regexes with some versions of gcc.
# NB pkg-config support for Boost is lacking; see https://svn.boost.org/trac/boost/ticket/1094
ifeq (,$(USING_BOOST))
BOOSTPREFIX =
BOOSTFLAGS =
BOOSTLIBS =
else
BOOSTPREFIX = /usr
ifeq (,$(wildcard $(BOOSTPREFIX)/include/boost/regex.h))
BOOSTPREFIX = /usr/local
ifeq (,$(wildcard $(BOOSTPREFIX)/include/boost/regex.h))
BOOSTPREFIX =
endif
endif

BOOSTFLAGS =
BOOSTLIBS =
ifneq (,$(BOOSTPREFIX))
BOOSTFLAGS := -DUSE_BOOST -I$(BOOSTPREFIX)/include
BOOSTLIBS := -L$(BOOSTPREFIX)/lib -lboost_regex
endif
endif

# install dir
PREFIX = /usr/local

# other flags
ifneq (,$(IS_DEBUG))
CPPFLAGS = -std=c++11 -g -DUSE_VECTOR_GUARDS $(GSLFLAGS) $(BOOSTFLAGS)
else
CPPFLAGS = -std=c++11 -g -O3 $(GSLFLAGS) $(BOOSTFLAGS)
endif
LIBFLAGS = -lstdc++ -lz $(GSLLIBS) $(BOOSTLIBS)

CPPFILES = $(wildcard src/*.cpp)
OBJFILES = $(subst src/,obj/,$(subst .cpp,.o,$(CPPFILES)))

# try clang++, fall back to g++
CPP = clang++
ifeq (, $(shell which $(CPP)))
CPP = g++
endif

# pwd
PWD = $(shell pwd)

# /bin/sh
SH = /bin/sh

# Targets

MAIN = historian

all: $(MAIN)

debug: $(MAIN)

$(MAIN): bin/$(MAIN)

install: bin/$(MAIN)
	cp $< $(PREFIX)/bin
	chmod a+x $(PREFIX)/bin/$(MAIN)

uninstall:
	rm $(PREFIX)/bin/$(MAIN)

clean:
	rm -rf bin/* obj/*

# Main build rules
bin/%: $(OBJFILES) obj/%.o
	@test -e bin || mkdir bin
	$(CPP) -o $@ obj/$*.o $(OBJFILES) $(LIBFLAGS)

obj/%.o: src/%.cpp
	@test -e obj || mkdir obj
	$(CPP) $(CPPFLAGS) -c -o $@ $<

obj/%.o: t/%.cpp
	@test -e obj || mkdir obj
	$(CPP) $(CPPFLAGS) -c -o $@ $<


# Tests

TEST = @perl/testexpect.pl
TEST4 = $(TEST) perl/roundfloats.pl 4
TEST10 = $(TEST) perl/roundfloats.pl 10

test: testregex testlogsumexp testseqio testnexus teststockholm testrateio testmatexp testmerge testseqprofile testforward testnullforward testbackward testnj testupgma testquickalign testtreeio testsubcount testnumsubcount testaligncount testsumprod testcountio testhist testcount testsum
# Skipped due to inconsistent platform-dependent behavior: testspan testhist-rndspan

testregex: bin/testregex
	$(TEST) bin/testregex /dev/null /dev/null

testlogsumexp: bin/testlogsumexp
	@bin/testlogsumexp -slow >data/logsumexp.txt 2> /dev/null
	$(TEST) bin/testlogsumexp -fast data/logsumexp.txt 2> /dev/null

testseqio: bin/testseqio
	$(TEST) bin/testseqio data/testaligncount.fa data/testaligncount.fa
	$(TEST) bin/testseqio data/gp120.fa data/gp120.fa

testnexus: bin/testnexus
	$(TEST) bin/testnexus data/testnexus.nex data/testnexus.nex

teststockholm: bin/teststockholm
	$(TEST) bin/teststockholm data/cbs.stock data/cbs.stock
	$(TEST) bin/teststockholm data/Lysine.stock data/Lysine.stock

testmatexp: bin/testmatexp
	$(TEST10) bin/testmatexp data/testrates.json 1 data/testrates.probs.json
	$(TEST10) bin/testmatexp -eigen data/testrates.json 1 data/testrates.probs.json

testrateio: bin/testrateio
	$(TEST4) bin/testrateio data/testrates.json data/testrates.out.json
	$(TEST4) bin/testrateio data/testrates.out.json data/testrates.out.json
	$(TEST4) bin/testrateio data/testrates.mix2.json data/testrates.mix2.out.json
	$(TEST4) bin/testrateio data/testrates.mix2.out.json data/testrates.mix2.out.json

testmerge: bin/testmerge
	$(TEST) bin/testmerge data/testmerge1.xy.fa data/testmerge1.xz.fa data/testmerge1.xyz.fa
	$(TEST) bin/testmerge data/testmerge1.xy.fa data/testmerge1.ayz.fa data/testmerge1.xyaz.fa
	$(TEST) bin/testmerge data/testmerge1.xz.fa data/testmerge1.ayz.fa data/testmerge1.xzay.fa
	$(TEST) bin/testmerge data/testmerge1.axyz.fa data/testmerge1.xz.fa data/testmerge1.axyz.fa
	$(TEST) bin/testmerge data/testmerge1.xy.fa data/testmerge1.xz.fa data/testmerge1-fail.ayz.fa data/empty 2> /dev/null
	$(TEST) bin/testmerge data/testmerge2.1.fa data/testmerge2.2.fa data/testmerge2.3.fa data/testmerge2.out.fa data/empty 2> /dev/null

testseqprofile: bin/testseqprofile
	$(TEST) bin/testseqprofile ACGT AAGCT data/testseqprofile.aagct.json

testforward: bin/testforward
	$(TEST) bin/testforward -all -matrix data/testforward.id100.len2.fa data/testforward.nosub.json 1 data/testforward.id100.len2.nosub.out
	@perl/testcumlp.pl data/testforward.id100.len2.nosub.out 51
	$(TEST) bin/testforward -hubs -best data/testforward.len2.fa data/testforward.nosub.json 1 data/testforward.len2.nosub.best.out
	$(TEST) bin/testforward -hubs -best data/testforward.len2.fa data/testforward.jukescantor.json 1 data/testforward.len2.jc.best.out
	$(TEST) bin/testforward -hubs -best data/testforward.len2-4.fa data/testforward.jukescantor.json .1 .01 data/testforward.len2-4.xdel.out
	$(TEST) bin/testforward -hubs -best data/testforward.len2-4.fa data/testforward.jukescantor.json .01 1 data/testforward.len2-4.yins.out
	$(TEST) bin/testforward -all 10 data/testforward.len2-4.fa data/testforward.jukescantor.json .1 data/testforward.len2-4.n10.all.out
	$(TEST) bin/testforward -hubs 10 data/testforward.len2-4.fa data/testforward.jukescantor.json .1 data/testforward.len2-4.n10.hubs.out

testnullforward: bin/testnullforward
	$(TEST) bin/testnullforward data/testforward.nosub.json 1 data/testnullforward.nosub.out

testbackward: bin/testbackward
	$(TEST) bin/testbackward data/testforward.len2.fa data/testforward.jukescantor.json 1 data/testbackward.len2.out
	$(TEST) bin/testbackward data/testforward.len2-4.fa data/testforward.jukescantor.json 1 data/testbackward.len2-4.out

testtreeio: bin/testtreeio
	$(TEST) bin/testtreeio data/PF16593.nhx data/PF16593.nhx
	$(TEST) bin/testtreeio data/testnj.out.nh data/testnj.out.nh
	$(TEST) bin/testtreeio data/PF16593.testspan.testnj.nh data/PF16593.testspan.testnj.nh
	$(TEST) bin/testtreeio data/testtreedupname.nh data/empty 2> /dev/null
	$(TEST) bin/testtreeio data/testtreenobranchlen.nh data/testtreenobranchlen.nh 2> /dev/null
	$(TEST) bin/testtreeio data/testreroot.nh C data/testreroot.c.nh

testspan: bin/testspan
	$(TEST) bin/testspan data/PF16593.fa data/testamino.json 1 data/PF16593.testspan.fa

testnj: bin/testnj
	$(TEST) bin/testnj data/testnj.jukescantor.json data/testnj.fa data/testnj.out.nh
	$(TEST) bin/testnj data/testamino.json data/PF16593.testspan.fa data/PF16593.testspan.testnj.nh

testupgma: bin/testupgma
	$(TEST) bin/testupgma data/testnj.jukescantor.json data/testnj.fa data/testupgma.out.nh
	$(TEST) bin/testupgma data/testamino.json data/PF16593.testspan.fa data/PF16593.testspan.testupgma.nh

testquickalign: bin/testquickalign
	$(TEST) bin/testquickalign data/PF16593.pair.fa data/testamino.json 1 data/testquickalign.out.fa

testsubcount: bin/testsubcount
	$(TEST) bin/testsubcount data/testrates.json A T 1 data/testsubcount1.json
	$(TEST) bin/testsubcount data/testforward.jukescantor.json A T 1 data/testsubcount2.json
	$(TEST) bin/testsubcount data/testforward.jukescantor.json A T 1 data/testsubcount3.json
	$(TEST) bin/testsubcount data/testrates.mix2.json A T 1 data/testsubcount.mix2.json

testnumsubcount: bin/testnumsubcount
	$(TEST) bin/testnumsubcount data/testforward.jukescantor.json A T A T .01 4 data/testnumsubcount1.out
	$(TEST) bin/testnumsubcount data/testforward.jukescantor.json A T A T 1 4 data/testnumsubcount2.out
	$(TEST) bin/testnumsubcount data/testforward.jukescantor.json A T C G 1 4 data/testnumsubcount3.out
	$(TEST) bin/testnumsubcount data/testrates.json A T A T 1 data/testnumsubcount4.out

testaligncount: bin/testaligncount
	$(TEST) bin/testaligncount -sub data/testnj.jukescantor.json data/testaligncount.fa data/testaligncount.nh data/testaligncount.out
	$(TEST) bin/testaligncount -eigen data/testnj.jukescantor.json data/testaligncount.fa data/testaligncount.nh data/testaligncount.out
	$(TEST) bin/testaligncount -sub data/testcount.jukescantor.json data/testaligncount2.fa data/testcount.nh data/testaligncount2.out.json

testsumprod: bin/testsumprod
	$(TEST) bin/testsumprod data/testnj.jukescantor.json data/testaligncount.fa data/testaligncount.nh data/testsumprod.out

testcountio: bin/testcountio
	$(TEST) bin/testcountio data/testcount.count.json data/testcount.count.json

testhist: bin/$(MAIN)
	$(TEST) bin/$(MAIN) recon -fast -norefine -output fasta -model data/testcount.jukescantor.json -guide data/testcount.fa -tree data/testcount.nh data/testcount.historian.fa
	$(TEST) bin/$(MAIN) recon -fast -norefine -output fasta -model data/testnj.jukescantor.json -nexus data/testnexus.nex data/testnexus.hist.fa
	$(TEST) bin/$(MAIN) recon -fast -norefine -output fasta -profsamples 100 -guide data/PF16593.testspan.fa -model data/testamino.json -tree data/PF16593.testspan.testnj.nh -band 10 data/PF16593.testspan.testnj.historian.fa
	$(TEST) bin/$(MAIN) recon -fast -norefine -output fasta -profsamples 100 -guide data/PF16593.testspan.fa -tree data/PF16593.testspan.testnj.nh -model data/testamino.json data/PF16593.testspan.testnj.historian.fa
	$(TEST) bin/$(MAIN) recon -fast -norefine -output fasta -profsamples 100 -guide data/PF16593.testspan.fa -model data/testamino.json -nj data/PF16593.testspan.testnj.historian.fa
	$(TEST) bin/$(MAIN) recon -fast -norefine -output fasta -profsamples 100 -seqs data/PF16593.fa -tree data/PF16593.nhx -model data/testamino.json -nj data/PF16593.historian.fa

testhist-rndspan:
	$(TEST) bin/$(MAIN) recon -fast -norefine -output fasta -profsamples 100 -rndspan data/PF16593.fa -model data/testamino.json -nj data/PF16593.testspan.testnj.historian.fa

testcount: bin/$(MAIN)
	$(TEST) bin/$(MAIN) count -fast -model data/testcount.jukescantor.json -recon data/testcount.fa -tree data/testcount.nh data/testcount.out.json
	$(TEST) bin/$(MAIN) count -fast -model data/testcount.jukescantor.json -tree data/testcount.nh -recon data/testcount.historian.fa data/testcount.count.json
	$(TEST) bin/$(MAIN) count -fast -model data/testrates.mix2.json -recon data/testcount.mix2.fa -tree data/testcount.mix2.nh data/testcount.mix2.count.json

testsum: bin/$(MAIN)
	$(TEST) bin/$(MAIN) sum data/testcount.out.json data/testcount.out.json data/testcount.sum.json

testgp120:
	bin/$(MAIN) recon -fast -norefine -guide data/gp120.guide.fa -tree data/gp120.tree.nh

testpost:
	bin/$(MAIN) post -fast -model data/testcount.jukescantor.json -guide data/testcount.fa -tree data/testcount.nh -v8

# Rules for building files in the repository
# For updating README.md
README.md: bin/$(MAIN)
	PATH=bin:$(PATH); $(MAIN) help | perl -pe 's/</&lt;/g;s/>/&gt;/g;' | perl -e 'open FILE,"<README.md";while(<FILE>){last if/<pre>/;print}close FILE;print"<pre><code>\n";while(<>){s/(default) of \S+ (uses at most \S+ of memory for DP matrix)/$$1 $$2/;print};print"</code></pre>\n"' >temp.md
	mv temp.md $@

# Rules for building release binaries
hide-shared:
	python/hide-shared-libs.py -d /usr/local --hide

restore-shared:
	python/hide-shared-libs.py -d /usr/local --restore
