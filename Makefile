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
	$(CPP) $(LIBFLAGS) -o $@ obj/$*.o $(OBJFILES)

obj/%.o: src/%.cpp
	@test -e obj || mkdir obj
	$(CPP) $(CPPFLAGS) -c -o $@ $<

obj/%.o: t/%.cpp
	@test -e obj || mkdir obj
	$(CPP) $(CPPFLAGS) -c -o $@ $<

# Tests

TEST = @perl/testexpect.pl

test: testlogsumexp testrateio testmatexp testmerge testseqprofile testforward testnullforward testnj testquickalign testspan testhist

testlogsumexp: bin/testlogsumexp
	@bin/testlogsumexp -slow >data/logsumexp.txt 2> /dev/null
	$(TEST) bin/testlogsumexp -fast data/logsumexp.txt 2> /dev/null

testmatexp: bin/testmatexp
	$(TEST) bin/testmatexp data/testrates.json 1 data/testrates.probs.json

testrateio: bin/testrateio
	$(TEST) bin/testrateio data/testrates.json data/testrates.out.json

testmerge: bin/testmerge
	$(TEST) bin/testmerge data/testmerge1.xy.fa data/testmerge1.xz.fa data/testmerge1.xyz.fa
	$(TEST) bin/testmerge data/testmerge1.xy.fa data/testmerge1.ayz.fa data/testmerge1.xyaz.fa
	$(TEST) bin/testmerge data/testmerge1.xz.fa data/testmerge1.ayz.fa data/testmerge1.xzay.fa
	$(TEST) bin/testmerge data/testmerge1.axyz.fa data/testmerge1.xz.fa data/testmerge1.axyz.fa
	$(TEST) bin/testmerge data/testmerge1.xy.fa data/testmerge1.xz.fa data/testmerge1-fail.ayz.fa data/empty 2> /dev/null

testseqprofile: bin/testseqprofile
	$(TEST) bin/testseqprofile ACGT AAGCT data/testseqprofile.aagct.json

testforward: bin/testforward
	$(TEST) bin/testforward -all -matrix data/testforward.id100.len2.fa data/testforward.nosub.json 1 data/testforward.id100.len2.nosub.out
	$(TEST) bin/testforward -absorbers -best data/testforward.len2.fa data/testforward.nosub.json 1 data/testforward.len2.nosub.best.out
	$(TEST) bin/testforward -absorbers -best data/testforward.len2.fa data/testforward.jukescantor.json 1 data/testforward.len2.jc.best.out
	$(TEST) bin/testforward -absorbers -best data/testforward.len2-4.fa data/testforward.jukescantor.json .1 .01 data/testforward.len2-4.xdel.out
	$(TEST) bin/testforward -absorbers -best data/testforward.len2-4.fa data/testforward.jukescantor.json .01 1 data/testforward.len2-4.yins.out
	$(TEST) bin/testforward -all 10 data/testforward.len2-4.fa data/testforward.jukescantor.json .1 data/testforward.len2-4.n10.all.out
	$(TEST) bin/testforward -absorbers 10 data/testforward.len2-4.fa data/testforward.jukescantor.json .1 data/testforward.len2-4.n10.abs.out
	$(TEST) bin/testforward -hubs 10 data/testforward.len2-4.fa data/testforward.jukescantor.json .1 data/testforward.len2-4.n10.hubs.out

testnullforward: bin/testnullforward
	$(TEST) bin/testnullforward data/testforward.nosub.json 1 data/testnullforward.nosub.out

testnj: bin/testnj
	$(TEST) bin/testnj data/testnj.jukescantor.json data/testnj.fa data/testnj.out.nh
	$(TEST) bin/testnj data/amino.json data/testspan.out.fa data/testnj.span.fa

testquickalign: bin/testquickalign
	$(TEST) bin/testquickalign data/PF16593.pair.fa data/amino.json 1 data/testquickalign.out.fa

testspan: bin/testspan
	$(TEST) bin/testspan data/PF16593.fa data/amino.json 1 data/testspan.out.fa

testhist: bin/historian
	$(TEST) bin/historian align -guide data/testspan.out.fa -model data/amino.json -tree data/testnj.span.fa -band 10 data/PF16593.testspan.testnj.historian.fa

testhist2: bin/historian
	bin/historian align -seqs data/PF16593.fa -tree data/PF16593.nhx -model data/amino.json -vv

# Rules for building files in the repository
# For updating README.md
README.md: bin/$(MAIN)
	PATH=bin:$(PATH); $(MAIN) help | perl -pe 's/</&lt;/g;s/>/&gt;/g;' | perl -e 'open FILE,"<README.md";while(<FILE>){last if/<pre>/;print}close FILE;print"<pre><code>\n";while(<>){print};print"</code></pre>\n"' >temp.md
	mv temp.md $@
