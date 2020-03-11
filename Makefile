# -*- makefile-gmake -*-

.SECONDARY:

# Pseudotargets that control compilation
USING_BOOST = $(findstring boost,$(MAKECMDGOALS))
IS_DEBUG = $(findstring debug,$(MAKECMDGOALS))
USING_EMSCRIPTEN = $(findstring emscripten,$(MAKECMDGOALS))

# If using emscripten, we need to compile gsl-js ourselves
ifneq (,$(USING_EMSCRIPTEN))
GSL_PREFIX = gsl-js
GSL_SOURCE = $(GSL_PREFIX)/gsl-js
GSL_LIB = $(GSL_PREFIX)/lib
GSL_FLAGS = -I$(GSL_SOURCE)
GSL_LIBS =
GSL_SUBDIRS = vector matrix utils linalg blas cblas block err min multimin permutation sys poly cdf complex eigen randist specfunc
GSL_OBJ_FILES = $(foreach dir,$(GSL_SUBDIRS),$(wildcard $(GSL_SOURCE)/$(dir)/*.o))
GSL_DEPS = $(GSL_LIB)
else
GSL_LIB =
GSL_DEPS =
GSL_OBJ_FILES =
# Try to figure out where GSL is
# autoconf would be better but we just need a quick hack for now :)
# Thanks to Torsten Seemann for gsl-config and pkg-config formulae
GSL_PREFIX = $(shell gsl-config --prefix)
ifeq (,$(wildcard $(GSL_PREFIX)/include/gsl/gsl_sf.h))
GSL_PREFIX = /usr
ifeq (,$(wildcard $(GSL_PREFIX)/include/gsl/gsl_sf.h))
GSL_PREFIX = /usr/local
endif
endif

GSL_FLAGS = $(shell pkg-config --cflags gsl)
ifeq (, $(GSL_FLAGS))
GSL_FLAGS = -I$(GSL_PREFIX)/include
endif

GSL_LIBS = $(shell pkg-config --libs gsl)
ifeq (, $(GSL_LIBS))
GSL_LIBS = -L$(GSL_PREFIX)/lib -lgsl -lgslcblas
endif
endif

# figure out whether to use Boost
# Boost is needed for regexes with some versions of gcc.
# NB pkg-config support for Boost is lacking; see https://svn.boost.org/trac/boost/ticket/1094
ifeq (,$(USING_BOOST))
BOOST_PREFIX =
BOOST_FLAGS =
BOOST_LIBS =
else
BOOST_PREFIX = /usr
ifeq (,$(wildcard $(BOOST_PREFIX)/include/boost/regex.h))
BOOST_PREFIX = /usr/local
ifeq (,$(wildcard $(BOOST_PREFIX)/include/boost/regex.h))
BOOST_PREFIX =
endif
endif

BOOST_FLAGS =
BOOST_LIBS =
ifneq (,$(BOOST_PREFIX))
BOOST_FLAGS := -DUSE_BOOST -I$(BOOST_PREFIX)/include
BOOST_LIBS := -L$(BOOST_PREFIX)/lib -lboost_regex
endif
endif

# install dir
PREFIX = /usr/local
INSTALL_BIN = $(PREFIX)/bin

# other flags
ALL_FLAGS = $(GSL_FLAGS) $(BOOST_FLAGS)
ALL_LIBS = $(GSL_LIBS) $(BOOST_LIBS)

ifneq (,$(IS_DEBUG))
CPP_FLAGS = -std=c++11 -g -DUSE_VECTOR_GUARDS -DDEBUG
else
ifneq (,$(IS_UNOPTIMIZED))
CPP_FLAGS = -std=c++11 -g
else
CPP_FLAGS = -std=c++11 -g -O3
endif
endif

ifneq (,$(USING_EMSCRIPTEN))
LD_FLAGS =
else
LD_FLAGS = -lstdc++ -lm
endif

CPP_FLAGS += $(ALL_FLAGS) -Isrc
LD_FLAGS += $(ALL_LIBS)

ifneq (,$(USING_EMSCRIPTEN))
EMCC_FLAGS = -s USE_ZLIB=1 -s EXTRA_EXPORTED_RUNTIME_METHODS="['FS', 'callMain']" -s ALLOW_MEMORY_GROWTH=1 -s EXIT_RUNTIME=1 --pre-js emcc/pre.js
CPP_FLAGS += $(EMCC_FLAGS)
LD_FLAGS += $(EMCC_FLAGS)
else
LD_FLAGS += -lz
endif

# files
CPP_FILES = $(wildcard src/*.cpp)
OBJ_FILES = $(subst src/,obj/,$(subst .cpp,.o,$(CPP_FILES)))

# C++ compiler: Emscripten, Clang, or GCC?
ifneq (,$(USING_EMSCRIPTEN))
CPP = emcc
else
# try clang++, fall back to g++
CPP = clang++
ifeq (, $(shell which $(CPP)))
CPP = g++
endif
endif

# pwd
PWD = $(shell pwd)

# /bin/sh
SH = /bin/sh

# Targets

MAIN = historian

ifneq (,$(USING_EMSCRIPTEN))
WRAP = node wasm/cmdwrap.js
MAINTARGET = wasm/historian.js
WRAPTARGET = $(WRAP) $(MAINTARGET)
TESTSUFFIX = .js
else
WRAP =
MAINTARGET = bin/$(MAIN)
WRAPTARGET = $(MAINTARGET)
TESTSUFFIX =
endif

all: $(MAIN)

debug: $(MAIN)

$(MAIN): bin/$(MAIN)

install: bin/$(MAIN)
	cp $< $(INSTALL_BIN)
	chmod a+x $(INSTALL_BIN)/$(MAIN)

uninstall:
	rm $(PREFIX)/bin/$(MAIN)

emscripten: $(MAINTARGET)

clean:
	rm -rf bin/* obj/*

# Main build rules
bin/% wasm/%.js: $(OBJ_FILES) obj/%.o $(GSL_DEPS)
	@test -e $(dir $@) || mkdir -p $(dir $@)
	$(CPP) $(LD_FLAGS) -o $@ obj/$*.o $(OBJ_FILES) $(GSL_OBJ_FILES)

obj/%.o: src/%.cpp $(GSL_DEPS)
	@test -e $(dir $@) || mkdir -p $(dir $@)
	$(CPP) $(CPP_FLAGS) -c -o $@ $<

obj/%.o: target/%.cpp $(GSL_DEPS)
	@test -e $(dir $@) || mkdir -p $(dir $@)
	$(CPP) $(CPP_FLAGS) -c -o $@ $<

bin/%: t/%.cpp $(OBJ_FILES)
	@test -e $(dir $@) || mkdir -p $(dir $@)
	$(CPP) $(CPP_FLAGS) -c -o obj/$*.o $<
	$(CPP) $(LD_FLAGS) -o $@$(TESTSUFFIX) obj/$*.o $(OBJ_FILES) $(GSL_OBJ_FILES)
	mv $@$(TESTSUFFIX) $@

# emscripten source files
# gsl-js
$(GSL_LIB):
	mkdir $(GSL_PREFIX)
	cd $(GSL_PREFIX); git clone https://github.com/GSL-for-JS/gsl-js.git
	cd $(GSL_SOURCE); emconfigure ./configure --prefix=$(abspath $(CURDIR)/$(GSL_PREFIX)); emmake make -k install

# Tests

TEST = @perl/testexpect.pl
WRAPTESTMAIN = $(TEST) $(WRAP) $(MAINTARGET)
WRAPTEST = $(TEST) $(WRAP)
WRAPTEST4 = $(TEST) perl/roundfloats.pl 4 $(WRAP)
WRAPTEST10 = $(TEST) perl/roundfloats.pl 10 $(WRAP)

test: testregex testlogsumexp testseqio testnexus teststockholm testrateio testmatexp testmerge testseqprofile testforward testnullforward testbackward testnj testupgma testquickalign testtreeio testsubcount testnumsubcount testaligncount testsumprod testcountio testhist testcount testsum
# Skipped due to inconsistent platform-dependent behavior: testspan testhist-rndspan

testregex: bin/testregex
	$(WRAPTEST) bin/testregex /dev/null /dev/null

testlogsumexp: bin/testlogsumexp
	@$(WRAP) bin/testlogsumexp -slow >data/logsumexp.txt 2> /dev/null
	$(WRAPTEST) bin/testlogsumexp -fast data/logsumexp.txt 2> /dev/null

testseqio: bin/testseqio
	$(WRAPTEST) bin/testseqio data/testaligncount.fa data/testaligncount.fa
	$(WRAPTEST) bin/testseqio data/gp120.fa data/gp120.fa

testnexus: bin/testnexus
	$(WRAPTEST) bin/testnexus data/testnexus.nex data/testnexus.nex

teststockholm: bin/teststockholm
	$(WRAPTEST) bin/teststockholm data/cbs.stock data/cbs.stock
	$(WRAPTEST) bin/teststockholm data/Lysine.stock data/Lysine.stock

testmatexp: bin/testmatexp
	$(WRAPTEST10) bin/testmatexp data/testrates.json 1 data/testrates.probs.json
	$(WRAPTEST10) bin/testmatexp -eigen data/testrates.json 1 data/testrates.probs.json

testrateio: bin/testrateio
	$(WRAPTEST4) bin/testrateio data/testrates.json data/testrates.out.json
	$(WRAPTEST4) bin/testrateio data/testrates.out.json data/testrates.out.json
	$(WRAPTEST4) bin/testrateio data/testrates.mix2.json data/testrates.mix2.out.json
	$(WRAPTEST4) bin/testrateio data/testrates.mix2.out.json data/testrates.mix2.out.json

testmerge: bin/testmerge
	$(WRAPTEST) bin/testmerge data/testmerge1.xy.fa data/testmerge1.xz.fa data/testmerge1.xyz.fa
	$(WRAPTEST) bin/testmerge data/testmerge1.xy.fa data/testmerge1.ayz.fa data/testmerge1.xyaz.fa
	$(WRAPTEST) bin/testmerge data/testmerge1.xz.fa data/testmerge1.ayz.fa data/testmerge1.xzay.fa
	$(WRAPTEST) bin/testmerge data/testmerge1.axyz.fa data/testmerge1.xz.fa data/testmerge1.axyz.fa
	$(WRAPTEST) bin/testmerge data/testmerge1.xy.fa data/testmerge1.xz.fa data/testmerge1-fail.ayz.fa data/empty 2> /dev/null
	$(WRAPTEST) bin/testmerge data/testmerge2.1.fa data/testmerge2.2.fa data/testmerge2.3.fa data/testmerge2.out.fa data/empty 2> /dev/null

testseqprofile: bin/testseqprofile
	$(WRAPTEST) bin/testseqprofile ACGT AAGCT data/testseqprofile.aagct.json

testforward: bin/testforward
	$(WRAPTEST) bin/testforward -all -matrix data/testforward.id100.len2.fa data/testforward.nosub.json 1 data/testforward.id100.len2.nosub.out
	@perl/testcumlp.pl data/testforward.id100.len2.nosub.out 51
	$(WRAPTEST) bin/testforward -hubs -best data/testforward.len2.fa data/testforward.nosub.json 1 data/testforward.len2.nosub.best.out
	$(WRAPTEST) bin/testforward -hubs -best data/testforward.len2.fa data/testforward.jukescantor.json 1 data/testforward.len2.jc.best.out
	$(WRAPTEST) bin/testforward -hubs -best data/testforward.len2-4.fa data/testforward.jukescantor.json .1 .01 data/testforward.len2-4.xdel.out
	$(WRAPTEST) bin/testforward -hubs -best data/testforward.len2-4.fa data/testforward.jukescantor.json .01 1 data/testforward.len2-4.yins.out
	$(WRAPTEST) bin/testforward -all 10 data/testforward.len2-4.fa data/testforward.jukescantor.json .1 data/testforward.len2-4.n10.all.out
	$(WRAPTEST) bin/testforward -hubs 10 data/testforward.len2-4.fa data/testforward.jukescantor.json .1 data/testforward.len2-4.n10.hubs.out

testnullforward: bin/testnullforward
	$(WRAPTEST) bin/testnullforward data/testforward.nosub.json 1 data/testnullforward.nosub.out

testbackward: bin/testbackward
	$(WRAPTEST) bin/testbackward data/testforward.len2.fa data/testforward.jukescantor.json 1 data/testbackward.len2.out
	$(WRAPTEST) bin/testbackward data/testforward.len2-4.fa data/testforward.jukescantor.json 1 data/testbackward.len2-4.out

testtreeio: bin/testtreeio
	$(WRAPTEST) bin/testtreeio data/PF16593.nhx data/PF16593.nhx
	$(WRAPTEST) bin/testtreeio data/testnj.out.nh data/testnj.out.nh
	$(WRAPTEST) bin/testtreeio data/PF16593.testspan.testnj.nh data/PF16593.testspan.testnj.nh
	$(WRAPTEST) bin/testtreeio data/testtreedupname.nh data/empty 2> /dev/null
	$(WRAPTEST) bin/testtreeio data/testtreenobranchlen.nh data/testtreenobranchlen.nh 2> /dev/null
	$(WRAPTEST) bin/testtreeio data/testreroot.nh C data/testreroot.c.nh

testspan: bin/testspan
	$(WRAPTEST) bin/testspan data/PF16593.fa data/testamino.json 1 data/PF16593.testspan.fa

testnj: bin/testnj
	$(WRAPTEST) bin/testnj data/testnj.jukescantor.json data/testnj.fa data/testnj.out.nh
	$(WRAPTEST) bin/testnj data/testamino.json data/PF16593.testspan.fa data/PF16593.testspan.testnj.nh

testupgma: bin/testupgma
	$(WRAPTEST) bin/testupgma data/testnj.jukescantor.json data/testnj.fa data/testupgma.out.nh
	$(WRAPTEST) bin/testupgma data/testamino.json data/PF16593.testspan.fa data/PF16593.testspan.testupgma.nh

testquickalign: bin/testquickalign
	$(WRAPTEST) bin/testquickalign data/PF16593.pair.fa data/testamino.json 1 data/testquickalign.out.fa

testsubcount: bin/testsubcount
	$(WRAPTEST) bin/testsubcount data/testrates.json A T 1 data/testsubcount1.json
	$(WRAPTEST) bin/testsubcount data/testforward.jukescantor.json A T 1 data/testsubcount2.json
	$(WRAPTEST) bin/testsubcount data/testforward.jukescantor.json A T 1 data/testsubcount3.json
	$(WRAPTEST) bin/testsubcount data/testrates.mix2.json A T 1 data/testsubcount.mix2.json

testnumsubcount: bin/testnumsubcount
	$(WRAPTEST) bin/testnumsubcount data/testforward.jukescantor.json A T A T .01 4 data/testnumsubcount1.out
	$(WRAPTEST) bin/testnumsubcount data/testforward.jukescantor.json A T A T 1 4 data/testnumsubcount2.out
	$(WRAPTEST) bin/testnumsubcount data/testforward.jukescantor.json A T C G 1 4 data/testnumsubcount3.out
	$(WRAPTEST) bin/testnumsubcount data/testrates.json A T A T 1 data/testnumsubcount4.out

testaligncount: bin/testaligncount
	$(WRAPTEST) bin/testaligncount -sub data/testnj.jukescantor.json data/testaligncount.fa data/testaligncount.nh data/testaligncount.out
	$(WRAPTEST) bin/testaligncount -eigen data/testnj.jukescantor.json data/testaligncount.fa data/testaligncount.nh data/testaligncount.out
	$(WRAPTEST) bin/testaligncount -sub data/testcount.jukescantor.json data/testaligncount2.fa data/testcount.nh data/testaligncount2.out.json

testsumprod: bin/testsumprod
	$(WRAPTEST) bin/testsumprod data/testnj.jukescantor.json data/testaligncount.fa data/testaligncount.nh data/testsumprod.out

testcountio: bin/testcountio
	$(WRAPTEST) bin/testcountio data/testcount.count.json data/testcount.count.json

testhist: $(MAINTARGET)
	$(WRAPTESTMAIN) recon -fast -norefine -output fasta -model data/testcount.jukescantor.json -guide data/testcount.fa -tree data/testcount.nh data/testcount.historian.fa
	$(WRAPTESTMAIN) recon -fast -norefine -output fasta -model data/testnj.jukescantor.json -nexus data/testnexus.nex data/testnexus.hist.fa
	$(WRAPTESTMAIN) recon -fast -norefine -output fasta -profsamples 100 -guide data/PF16593.testspan.fa -model data/testamino.json -tree data/PF16593.testspan.testnj.nh -band 10 data/PF16593.testspan.testnj.historian.fa
	$(WRAPTESTMAIN) recon -fast -norefine -output fasta -profsamples 100 -guide data/PF16593.testspan.fa -tree data/PF16593.testspan.testnj.nh -model data/testamino.json data/PF16593.testspan.testnj.historian.fa
	$(WRAPTESTMAIN) recon -fast -norefine -output fasta -profsamples 100 -guide data/PF16593.testspan.fa -model data/testamino.json -nj data/PF16593.testspan.testnj.historian.fa
	$(WRAPTESTMAIN) recon -fast -norefine -output fasta -profsamples 100 -seqs data/PF16593.fa -tree data/PF16593.nhx -model data/testamino.json -nj data/PF16593.historian.fa

testhist-rndspan:
	$(WRAPTESTMAIN) recon -fast -norefine -output fasta -profsamples 100 -rndspan data/PF16593.fa -model data/testamino.json -nj data/PF16593.testspan.testnj.historian.fa

testcount: $(MAINTARGET)
	$(WRAPTESTMAIN) count -fast -model data/testcount.jukescantor.json -recon data/testcount.fa -tree data/testcount.nh data/testcount.out.json
	$(WRAPTESTMAIN) count -fast -model data/testcount.jukescantor.json -tree data/testcount.nh -recon data/testcount.historian.fa data/testcount.count.json
	$(WRAPTESTMAIN) count -fast -model data/testrates.mix2.json -recon data/testcount.mix2.fa -tree data/testcount.mix2.nh data/testcount.mix2.count.json

testsum: $(MAINTARGET)
	$(WRAPTESTMAIN) sum data/testcount.out.json data/testcount.out.json data/testcount.sum.json

testgp120:
	$(MAINTARGET) recon -fast -norefine -guide data/gp120.guide.fa -tree data/gp120.tree.nh

testpost:
	$(MAINTARGET) post -fast -model data/testcount.jukescantor.json -guide data/testcount.fa -tree data/testcount.nh -v8

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
