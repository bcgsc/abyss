#!/usr/bin/make -Rrf

# Test abyss-connectpairs

SHELL=/bin/bash -o pipefail

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
# shared options to connectpairs
cp_opts=-j$j -v -v -k$k
# user options to connectpairs
CP_OPTS?=

#------------------------------------------------------------
# phony targets
#------------------------------------------------------------

tests=run_test \
	save_and_load_test \
	union_test \
	interleaved_files_test

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
	/usr/bin/time -v $(connectpairs) $(cp_opts) -b$b -o $(tmpdir)/e$* $(CP_OPTS) $^

$(tmpdir)/e%_l2.bloom: $(tmpdir) $(tmpdir)/e%_1.fq $(tmpdir)/e%_2.fq
	$(bloom) build -v -k$k -l2 -b$b $@ $(filter-out $<, $^)

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

save_and_load_test: $(tmpdir)/e$e_l2.bloom \
	$(tmpdir)/e$e_merged.fa \
	$(tmpdir)/e$e_reads_1.fq \
	$(tmpdir)/e$e_reads_2.fq
	/usr/bin/time -v $(connectpairs) $(cp_opts) -o $(tmpdir)/e$e_loaded \
		-i $(tmpdir)/e$e_l2.bloom $(CP_OPTS) $(tmpdir)/e$e_1.fq $(tmpdir)/e$e_2.fq
	diff $(tmpdir)/e$e_merged.fa $(tmpdir)/e$e_loaded_merged.fa
	diff $(tmpdir)/e$e_reads_1.fq $(tmpdir)/e$e_loaded_reads_1.fq
	diff $(tmpdir)/e$e_reads_2.fq $(tmpdir)/e$e_loaded_reads_2.fq
	@echo '------------------'
	@echo '$@: PASSED'
	@echo '------------------'

#------------------------------------------------------------
# interleaved_files_test
#------------------------------------------------------------

HALF_FASTQ_LINES:=$(shell echo '$N * 2 * 4 / 2' | bc)

interleaved_files_test: $(tmpdir)/e$e_l2.bloom \
		$(tmpdir)/e$e_merged.fa \
		$(tmpdir)/e$e_interleaved_a.fq \
		$(tmpdir)/e$e_interleaved_b.fq
	/usr/bin/time -v $(connectpairs) $(cp_opts) -I -b$b \
		-i $(tmpdir)/e$e_l2.bloom -o $(tmpdir)/e$e_interleaved \
		$(CP_OPTS) \
		$(tmpdir)/e$e_interleaved_a.fq \
		$(tmpdir)/e$e_interleaved_b.fq
	diff $(tmpdir)/e$e_merged.fa $(tmpdir)/e$e_interleaved_merged.fa
	diff $(tmpdir)/e$e_reads_1.fq $(tmpdir)/e$e_interleaved_reads_1.fq
	diff $(tmpdir)/e$e_reads_2.fq $(tmpdir)/e$e_interleaved_reads_2.fq
	@echo '------------------'
	@echo '$@: PASSED'
	@echo '------------------'

$(tmpdir)/e%_interleaved.fq: $(tmpdir)/e%_1.fq $(tmpdir)/e%_2.fq
	paste -d'\n' <(cat $(tmpdir)/e$e_1.fq | paste - - - -) \
		<(cat $(tmpdir)/e$e_2.fq | paste - - - -) | \
		tr '\t' '\n' > $(tmpdir)/e$e_interleaved.fq

#split -a 1 --additional-suffix .fq -n l/2 -d $< $(tmpdir)/e$*_interleaved_
$(tmpdir)/e%_interleaved_a.fq $(tmpdir)/e%_interleaved_b.fq: \
	$(tmpdir)/e%_interleaved.fq
	head -n $(HALF_FASTQ_LINES) $< > $(tmpdir)/e$*_interleaved_a.fq
	tail -n $(HALF_FASTQ_LINES) $< > $(tmpdir)/e$*_interleaved_b.fq

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
