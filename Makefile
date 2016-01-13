.SECONDARY:

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
# Boost is optional -- it's only needed for regexes with gcc
# NB pkg-config support for Boost is lacking; see https://svn.boost.org/trac/boost/ticket/1094
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

# install dir
PREFIX = /usr/local

# other flags
CPPFLAGS = -DUSE_VECTOR_GUARDS -std=c++11 -g $(GSLFLAGS) $(BOOSTFLAGS)
LIBFLAGS = -lstdc++ -lz $(GSLLIBS) $(BOOSTLIBS)

CPPFILES = $(wildcard src/*.cpp)

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

MAIN = propalor

all: $(MAIN)

$(MAIN): bin/$(MAIN)

install: bin/$(MAIN)
	cp $< $(PREFIX)/bin
	chmod a+x $(PREFIX)/bin/$(MAIN)

uninstall:
	rm $(PREFIX)/bin/$(MAIN)

clean:
	rm -rf bin/*

# Main build rule
bin/%: $(CPPFILES) t/%.cpp
	test -e bin || mkdir bin
	$(CPP) $(CPPFLAGS) $(LIBFLAGS) -o $@ t/$*.cpp $(CPPFILES)

# Tests

test: testlogsumexp

testlogsumexp: bin/testlogsumexp
	bin/testlogsumexp -slow >data/logsumexp.txt
	perl/testexpect.pl bin/testlogsumexp -fast data/logsumexp.txt


# quaff tests (integration tests)
# Tests of basic commands
quaff-tests: testquaffcountself testquaffalignself testquaffoverlapself testquaffcountself-remote testquaffalignself-remote testquaffoverlapself-remote testquaffcountself-qsub testquaffalignself-qsub testquaffoverlapself-qsub

testquaffcountself: bin/quaff
	perl/testexpect.pl bin/quaff count data/c8f30.fastq.gz data/c8f30.fastq.gz -kmatchmb 10 -fwdstrand data/c8f30-self-counts.json

testquaffalignself: bin/quaff
	perl/testexpect.pl bin/quaff align data/c8f30.fastq.gz data/c8f30.fastq.gz -kmatchmb 10 -fwdstrand data/c8f30-self-align.json

testquaffoverlapself: bin/quaff data/copy-of-c8f30.fastq
	perl/testexpect.pl bin/quaff overlap data/c8f30.fastq.gz data/copy-of-c8f30.fastq -kmatchmb 10 -fwdstrand data/c8f30-self-overlap.json

data/copy-of-c8f30.fastq: data/c8f30.fastq.gz
	gzcat $< | perl -pe s/channel/copy/ >$@

# Tests of the -remote option (parallelization over sockets)
testquaffcountself-remote: bin/quaff
	perl/testexpect.pl bin/quaff count $(PWD)/data/c8f30.fastq.gz $(PWD)/data/c8f30.fastq.gz -kmatchmb 10 -fwdstrand -remotepath $(PWD)/bin/quaff -remote localhost:$(RANDOMPORT) data/c8f30-self-counts.json

testquaffalignself-remote: bin/quaff
	perl/testexpect.pl bin/quaff align $(PWD)/data/c8f30.fastq.gz $(PWD)/data/c8f30.fastq.gz -kmatchmb 10 -fwdstrand -remotepath $(PWD)/bin/quaff -remote localhost:$(RANDOMPORT) data/c8f30-self-align.json

testquaffoverlapself-remote: bin/quaff data/copy-of-c8f30.fastq
	perl/testexpect.pl bin/quaff overlap $(PWD)/data/c8f30.fastq.gz $(PWD)/data/copy-of-c8f30.fastq -kmatchmb 10 -fwdstrand -remotepath $(PWD)/bin/quaff -remote localhost:$(RANDOMPORT) data/c8f30-self-overlap.json

# Test of -qsub* options (parallelization via shell scripts and NFS coordination; a.k.a. Sun Grid Engine, PBS and the like)
testquaffcountself-qsub: bin/quaff
	perl/testexpect.pl bin/quaff count $(PWD)/data/c8f30.fastq.gz $(PWD)/data/c8f30.fastq.gz -kmatchmb 10 -fwdstrand -remotepath $(PWD)/bin/quaff -qsubjobs 1 -qsubpath $(SH) data/c8f30-self-counts.json

testquaffalignself-qsub: bin/quaff
	perl/testexpect.pl bin/quaff align $(PWD)/data/c8f30.fastq.gz $(PWD)/data/c8f30.fastq.gz -kmatchmb 10 -fwdstrand -remotepath $(PWD)/bin/quaff -qsubjobs 1 -qsubpath $(SH) data/c8f30-self-align.json

testquaffoverlapself-qsub: bin/quaff data/copy-of-c8f30.fastq
	perl/testexpect.pl bin/quaff overlap $(PWD)/data/c8f30.fastq.gz $(PWD)/data/copy-of-c8f30.fastq -kmatchmb 10 -fwdstrand -remotepath $(PWD)/bin/quaff -qsubjobs 1 -qsubpath $(SH) data/c8f30-self-overlap.json




# testaws should be run separately (it needs secret keys, could cost money, etc)
testaws: bin/testaws
	echo Run like this: 'bin/testaws keyPairName'


# Rules for building files in the repository
# For updating README.md
README.md: bin/quaff
	PATH=bin:$(PATH); quaff help | perl -pe 's/</&lt;/g;s/>/&gt;/g;' | perl -e 'open FILE,"<README.md";while(<FILE>){last if/<pre>/;print}close FILE;print"<pre><code>\n";while(<>){print};print"</code></pre>\n"' >temp.md
	mv temp.md $@

# For updating default params
src/defaultparams.cpp: data/defaultparams.json
	perl -e 'open S,"<".shift();while(<S>){print;last if/defaultQuaffParamsText =/}close S;open A,"<".shift();$$q=chr(34);while(<A>){chomp;s/$$q/\\$$q/g;print chr(34),$$_,"\\n",chr(34),"\n"}print";\n"' $@ $< >temp.cpp
	mv temp.cpp $@
