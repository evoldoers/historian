
#MAKEDIR := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
MAKEDIR := .
#HISTDIR := $(MAKEDIR)/../..
HISTDIR := $(HOME)/historian

HISTPROG = $(HISTDIR)/bin/historian
HISTARGS = -v3 -output fasta -noancs
HISTORIAN = $(HISTPROG) $(HISTARGS)

QSCORE = qscore/qscore

PROT1 = $(HISTDIR)/model/prot1.json
PROT3 = $(HISTDIR)/model/prot3.json
PROT4 = $(HISTDIR)/model/prot4.json
PROT14 = $(HISTDIR)/model/prot1-4.json

ALL_FILES = $(notdir $(wildcard in/*))

# top level
METHOD ?= nj
ALL = $(addprefix historian/$(METHOD)/qscore/,$(ALL_FILES))

all: $(ALL)

DOT = $(addprefix historian/dot/,$(ALL_FILES))
dot:
	make $(DOT)

# Download benchmark & qscore
bench1.0:
	curl -O http://www.drive5.com/bench/bench.tar.gz
	tar xvzf bench.tar.gz
	cp Makefile $@

$(QSCORE):
	test -e qscore || mkdir qscore
	cd qscore; curl -O http://www.drive5.com/qscore/qscore_src.tar.gz
	cd qscore; tar xvzf qscore_src.tar.gz; make

# convenient pseudotarget
%: historian/$(METHOD)/qscore/%
	 make METHOD=$(METHOD) historian/$(METHOD)/qscore/$@

# directories
historian/$(METHOD)/align historian/$(METHOD)/qscore historian/$(METHOD)/guide historian/$(METHOD)/params historian/dot:
	mkdir -p $@

# qscore
historian/$(METHOD)/qscore/%: historian/$(METHOD)/align/% historian/$(METHOD)/qscore
	$(QSCORE) -seqdiffwarn -ignoremissingseqs -ignoretestcase -test historian/$(METHOD)/align/$* -ref ref/$* -truncname -quiet >$@

# nj
historian/nj/align/%: in/% historian/nj/align
	$(HISTORIAN) -fast in/$* -model $(PROT1) -refine -nj >$@

# fit
historian/fit/align/%: in/% historian/fit/align historian/fit/guide historian/fit/params
	$(HISTPROG) -fast in/$* -nj -v3 -output fasta -noancs >historian/fit/guide/$*
	$(HISTPROG) fit -fast -v3 -model $(PROT1) -fixsubrates -guide historian/fit/guide/$* >historian/fit/params/$*
	$(HISTORIAN) -fast -model historian/fit/params/$* -guide historian/fit/guide/$* -norefine >historian/fit/align/$*

# mcmc
historian/mcmc/align/%: in/% historian/mcmc/align
	$(HISTORIAN) -fast in/$* -model $(PROT1) -refine -mcmc >$@

# prot3
historian/prot3/align/%: in/% historian/prot3/align
	$(HISTORIAN) -fast in/$* -refine -nj -model $(PROT3) >$@

# prot4
historian/prot4/align/%: in/% historian/prot4/align
	$(HISTORIAN) -fast in/$* -refine -nj -model $(PROT4) >$@

# prot14
historian/prot1-4/align/%: in/% historian/prot1-4/align
	$(HISTORIAN) -fast in/$* -refine -nj -model $(PROT14) >$@

# fit4
historian/fit4/align/%: in/% historian/fit4/align historian/fit4/guide historian/fit4/params
	$(HISTPROG) -fast in/$* -nj -v3 >historian/fit4/guide/$*
	$(HISTPROG) -fast fit -v3 -model $(PROT4) -fixsubrates -stockholm historian/fit4/guide/$* >historian/fit4/params/$*
	$(HISTORIAN) -fast -model historian/fit4/params/$* -stockholm historian/fit4/guide/$* >historian/fit4/align/$*

# wag
historian/wag/align/%: in/% historian/wag/align
	$(HISTORIAN) -fast in/$* -refine -nj -preset wag >$@

# lg
historian/lg/align/%: in/% historian/lg/align
	$(HISTORIAN) -fast in/$* -refine -nj -preset lg >$@

# lg, slow
historian/lgslow/align/%: in/% historian/lgslow/align
	$(HISTORIAN) in/$* -refine -nj -preset lg >$@

# lg, faster
historian/lgfaster/align/%: in/% historian/lgfaster/align
	$(HISTORIAN) -faster in/$* -nj -preset lg >$@

# wag4
historian/wag4/align/%: in/% historian/wag4/align
	$(HISTORIAN) -fast in/$* -refine -nj -preset wag -gamma 4 >$@

# lg4
historian/lg4/align/%: in/% historian/lg4/align
	$(HISTORIAN) -fast in/$* -refine -nj -preset lg -gamma 4 >$@

# lg4, slow
historian/lg4slow/align/%: in/% historian/lg4slow/align
	$(HISTORIAN) in/$* -refine -nj -preset lg -gamma 4 >$@

# dot
historian/dot/%: in/% historian/dot
	$(HISTORIAN) -fast in/$* -savedot $@

# table
GETMEANS = perl -e 'use List::Util qw(sum);my(@q,@tc);while(<>){if(/Q=([\d\.]+);TC=([\d\.]+)/){push@q,$$1;push@tc,$$2}}print" & ",sum(@q)/@q," & ",sum(@tc)/@tc'
table.tex:
	@echo "Historian"
	@cat bali3/historian/nj/qscore/* | $(GETMEANS)
#	@cat prefab4/historian/nj/qscore/* | $(GETMEANS)
	@echo " \\\\\\\\"
	@echo "Historian (WAG)"
	@cat bali3/historian/wag/qscore/* | $(GETMEANS)
#	@cat prefab4/historian/wag/qscore/* | $(GETMEANS)
	@echo " \\\\\\\\"
	@echo "Historian (LG)"
	@cat bali3/historian/lg/qscore/* | $(GETMEANS)
#	@cat prefab4/historian/lg/qscore/* | $(GETMEANS)
	@echo " \\\\\\\\"
	@echo "Historian (LG, slow)"
	@cat bali3/historian/lgslow/qscore/* | $(GETMEANS)
	@echo " \\\\\\\\"
	@echo "CLUSTALW"
	@cat bali3/qscore/clw | $(GETMEANS)
#	@cat prefab4/qscore/clw | $(GETMEANS)
	@echo " \\\\\\\\"
	@echo "PRANK +F, CLUSTALW tree"
	@cat bali3/qscore/prankfclwtree | $(GETMEANS)
#	@cat prefab4/qscore/prankfclwtree | $(GETMEANS)
	@echo " \\\\\\\\"
	@echo "PRANK +F"
	@cat bali3/prank/qscore/* | $(GETMEANS)
	@echo " \\\\\\\\"
	@echo "MUSCLE v3.8.31"
	@cat bali3/muscle/qscore/* | $(GETMEANS)
	@echo " \\\\\\\\"
	@echo "MUSCLE v4.0"
	@cat bali3/qscore/mus4 | $(GETMEANS)
#	@cat prefab4/qscore/mus4 | $(GETMEANS)
	@echo " \\\\\\\\"
	@echo "MAFFT v6.603"
	@cat bali3/qscore/mafft6.603einsi | $(GETMEANS)
#	@cat prefab4/qscore/mafft6.603einsi | $(GETMEANS)
	@echo " \\\\\\\\"
	@echo "ProbCons 1.12"
	@cat bali3/qscore/pc1.12 | $(GETMEANS)
#	@cat prefab4/qscore/pc1.12 | $(GETMEANS)
	@echo " \\\\\\\\"



# Muscle
MUSCLE := muscle

ALLMUSCLE = $(addprefix muscle/qscore/,$(ALL_FILES))

allmuscle: $(ALLMUSCLE)

muscle/qscore/%: in/%
	$(MUSCLE) -in in/$* >muscle/align/$*
	$(QSCORE) -seqdiffwarn -ignoremissingseqs -ignoretestcase -test muscle/align/$* -ref ref/$* -truncname -quiet >muscle/qscore/$*


# Prank

PRANK := prank

ALLPRANK = $(addprefix prank/qscore/,$(ALL_FILES))

allprank: $(ALLPRANK)

prank/qscore/%: in/%
	$(PRANK) -d=in/$* -noxml -realbranches -writeanc -o=prank/output/$* +F
	cp prank/output/$*.1.fas prank/align/$*
	$(QSCORE) -seqdiffwarn -ignoremissingseqs -ignoretestcase -test prank/align/$* -ref ref/$* -truncname -quiet >prank/qscore/$*

# don't delete intermediate files
.SECONDARY:
