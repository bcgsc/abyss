#!/usr/bin/make -Rrf

# Test abyss-connectpairs

#------------------------------------------------------------
# testing params
#------------------------------------------------------------

# threads (leave as 1 for consistent results)
j?=1
# counting bloom filter size
b?=100M
# kmer size
k?=20
# number of synthetic read pairs
N?=5000
# error rate of synthetic reads
e?=0.005
# path to abyss-connectpairs binary
connectpairs?=abyss-connectpairs
# path to abyss-bloom binary
bloom?=abyss-bloom
# temp dir for test outputs
tmpdir=tmp

#------------------------------------------------------------
# phony targets
#------------------------------------------------------------

tests=run_test \
	save_and_load_test \
	union_test

.PHONY: all $(tests)
.DELETE_ON_ERROR:
.SECONDARY:

#------------------------------------------------------------
# top level rules
#------------------------------------------------------------

all: $(tests)

clean:
	rm -rf $(tmpdir)

#------------------------------------------------------------
# common rules
#------------------------------------------------------------

$(tmpdir):
	mkdir -p $@

$(tmpdir)/test_reference.fa: $(tmpdir)
	curl https://raw.github.com/dzerbino/velvet/master/data/test_reference.fa \
		|abyss-tofastq --fasta >$@

$(tmpdir)/e%_1.fq $(tmpdir)/e%_2.fq: $(tmpdir)/test_reference.fa
	wgsim -S 0 -e $* -N $N -r 0 -R 0 $< $(tmpdir)/e$*_1.fq $(tmpdir)/e$*_2.fq

$(tmpdir)/e%_merged.fa $(tmpdir)/e%_reads_1.fq $(tmpdir)/e%_reads_2.fq: $(tmpdir)/e%_1.fq $(tmpdir)/e%_2.fq
	/usr/bin/time -v $(connectpairs) -j$j -v -v -b$b -k$k -o $(tmpdir)/e$* $^

#------------------------------------------------------------
# run_test
#------------------------------------------------------------

run_test: $(tmpdir) $(tmpdir)/e$e_merged.fa
	@echo '------------------'
	@echo '$@: PASSED'
	@echo '------------------'

#------------------------------------------------------------
# save_and_load_test
#------------------------------------------------------------

save_and_load_test: $(tmpdir) \
	$(tmpdir)/e$e_merged.fa \
	$(tmpdir)/e$e_reads_1.fq \
	$(tmpdir)/e$e_reads_2.fq
	$(bloom) build -v -k$k -l2 -b$b $(tmpdir)/e$e_l2.bloom \
		$(tmpdir)/e$e_1.fq $(tmpdir)/e$e_2.fq
	/usr/bin/time -v $(connectpairs) -j$j -v -v -k$k -o $(tmpdir)/e$e_loaded \
		-i $(tmpdir)/e$e_l2.bloom  $(tmpdir)/e$e_1.fq $(tmpdir)/e$e_2.fq
	diff $(tmpdir)/e$e_merged.fa $(tmpdir)/e$e_loaded_merged.fa
	diff $(tmpdir)/e$e_reads_1.fq $(tmpdir)/e$e_loaded_reads_1.fq
	diff $(tmpdir)/e$e_reads_2.fq $(tmpdir)/e$e_loaded_reads_2.fq
	@echo '------------------'
	@echo '$@: PASSED'
	@echo '------------------'

#------------------------------------------------------------
# union_test
#------------------------------------------------------------

union_test: $(tmpdir) $(tmpdir)/e$e_1.fq $(tmpdir)/e$e_2.fq
	$(bloom) build -v -k$k -b$b $(tmpdir)/e$e.bloom $(tmpdir)/e$e_1.fq $(tmpdir)/e$e_2.fq
	$(bloom) build -v -k$k -b$b $(tmpdir)/e$e_1.bloom $(tmpdir)/e$e_1.fq
	$(bloom) build -v -k$k -b$b $(tmpdir)/e$e_2.bloom $(tmpdir)/e$e_2.fq
	$(bloom) union -v -k$k $(tmpdir)/e$e_union.bloom $(tmpdir)/e$e_1.bloom $(tmpdir)/e$e_2.bloom
	cmp $(tmpdir)/e$e.bloom $(tmpdir)/e$e_union.bloom
	@echo '------------------'
	@echo '$@: PASSED'
	@echo '------------------'
