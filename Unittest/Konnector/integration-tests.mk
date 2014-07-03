#!/usr/bin/make -Rrf

# Test konnector

SHELL=/bin/bash -o pipefail

#------------------------------------------------------------
# testing params
#------------------------------------------------------------

# threads (leave as 1 for consistent results)
j?=1
# counting bloom filter size
# For the parallel_load_* tests, this number
# must be in bytes (no "M" or "G" suffix) and divisible by 2.
b?=100000000
# kmer size
k?=20
# number of synthetic read pairs
# For the parallel_load_* tests, this number must be
# divisible by 3.
N?=6000
# error rate of synthetic reads
e?=0.005
# path to konnector binary
konnector?=konnector
# path to abyss-bloom binary
bloom?=abyss-bloom
# path to abyss-bloom-dist.mk makefile
bloom_dist=abyss-bloom-dist.mk
# temp dir for test outputs
tmpdir=tmp
# shared options to konnector
k_opts=-j$j -v -v -k$k
# user options to konnector
K_OPTS?=

#------------------------------------------------------------
# global vars
#------------------------------------------------------------

b_div_2:=$(shell echo $b / 2 | bc)

#------------------------------------------------------------
# phony targets
#------------------------------------------------------------

tests=run_test \
	save_and_load_test \
	union_test \
	interleaved_files_test \
	window_test \
	parallel_load_2_files_test \
	parallel_load_3_files_test \
	abyss_bloom_dist_1_file_test \
	abyss_bloom_dist_2_files_test \
	abyss_bloom_dist_3_files_test \
	abyss_bloom_illegal_chars_test

.PHONY: all $(tests)
.DELETE_ON_ERROR:
.SECONDARY:

#------------------------------------------------------------
# top level rules
#------------------------------------------------------------

.PHONY: all clean tmpdir

all: $(tests)

clean:
	rm -f $(tmpdir)/*
	rmdir $(tmpdir) || true

#------------------------------------------------------------
# common rules
#------------------------------------------------------------

$(tmpdir):
	mkdir -p $(tmpdir)

$(tmpdir)/test_reference.fa: | $(tmpdir)
	curl -L https://raw.github.com/dzerbino/velvet/master/data/test_reference.fa \
		|abyss-tofastq --fasta >$@

$(tmpdir)/e%_1.fq $(tmpdir)/e%_2.fq: $(tmpdir)/test_reference.fa
	wgsim -S 0 -e $* -N $N -r 0 -R 0 $< $(tmpdir)/e$*_1.fq $(tmpdir)/e$*_2.fq

$(tmpdir)/e%_merged.fa $(tmpdir)/e%_reads_1.fq $(tmpdir)/e%_reads_2.fq: $(tmpdir)/e%_1.fq $(tmpdir)/e%_2.fq
	/usr/bin/time -v $(konnector) $(k_opts) -b$b -o $(tmpdir)/e$* $(K_OPTS) $^

$(tmpdir)/e%_l2.bloom: $(tmpdir) $(tmpdir)/e%_1.fq $(tmpdir)/e%_2.fq
	$(bloom) build -v -k$k -l2 -b$b $@ $(filter-out $<, $^)

$(tmpdir)/e%_interleaved.fq: $(tmpdir)/e%_1.fq $(tmpdir)/e%_2.fq
	paste -d'\n' <(cat $(tmpdir)/e$e_1.fq | paste - - - -) \
		<(cat $(tmpdir)/e$e_2.fq | paste - - - -) | \
		tr '\t' '\n' > $(tmpdir)/e$e_interleaved.fq

FASTQ_CHUNKS:=3
FASTQ_CHUNK_SIZE:=$(shell echo '$N * 2 * 4 / $(FASTQ_CHUNKS)' | bc)

$(tmpdir)/e%_reads_1of3.fq \
	$(tmpdir)/e%_reads_2of3.fq \
	$(tmpdir)/e%_reads_3of3.fq: $(tmpdir)/e%_interleaved.fq
	awk '{ print > "$(tmpdir)/e$e_reads_" \
		int((NR-1)/$(FASTQ_CHUNK_SIZE))+1 \
		"of"$(FASTQ_CHUNKS)".fq"}' $<

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
	/usr/bin/time -v $(konnector) $(k_opts) -o $(tmpdir)/e$e_loaded \
		-i $(tmpdir)/e$e_l2.bloom $(K_OPTS) $(tmpdir)/e$e_1.fq $(tmpdir)/e$e_2.fq
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
	/usr/bin/time -v $(konnector) $(k_opts) -I -b$b \
		-i $(tmpdir)/e$e_l2.bloom -o $(tmpdir)/e$e_interleaved \
		$(K_OPTS) \
		$(tmpdir)/e$e_interleaved_a.fq \
		$(tmpdir)/e$e_interleaved_b.fq
	diff $(tmpdir)/e$e_merged.fa $(tmpdir)/e$e_interleaved_merged.fa
	diff $(tmpdir)/e$e_reads_1.fq $(tmpdir)/e$e_interleaved_reads_1.fq
	diff $(tmpdir)/e$e_reads_2.fq $(tmpdir)/e$e_interleaved_reads_2.fq
	@echo '------------------'
	@echo '$@: PASSED'
	@echo '------------------'

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

#------------------------------------------------------------
# intersect_test
#------------------------------------------------------------

intersect_test: $(tmpdir) $(tmpdir)/e$e_1.fq $(tmpdir)/e$e_2.fq
	$(bloom) build -v -k$k -b$b $(tmpdir)/e$e.bloom $(tmpdir)/e$e_1.fq $(tmpdir)/e$e_2.fq
	$(bloom) build -v -k$k -b$b $(tmpdir)/e$e_1.bloom $(tmpdir)/e$e_1.fq
	$(bloom) intersect -v -k$k $(tmpdir)/e$e_intersect.bloom \
		$(tmpdir)/e$e.bloom $(tmpdir)/e$e_1.bloom
	cmp $(tmpdir)/e$e_intersect.bloom $(tmpdir)/e$e_1.bloom
	@echo '------------------'
	@echo '$@: PASSED'
	@echo '------------------'

#------------------------------------------------------------
# window_test
#------------------------------------------------------------

window_test: $(tmpdir)/e$e_l2.bloom $(tmpdir)/e$e_1.fq $(tmpdir)/e$e_2.fq
	$(bloom) build -v -k$k -l2 -w1/2 -b$b $(tmpdir)/e$e_l2_window1.bloom \
		$(tmpdir)/e$e_1.fq $(tmpdir)/e$e_2.fq
	$(bloom) build -v -k$k -l2 -w2/2 -b$b $(tmpdir)/e$e_l2_window2.bloom \
		$(tmpdir)/e$e_1.fq $(tmpdir)/e$e_2.fq
	$(bloom) union -v -k$k $(tmpdir)/e$e_l2_concat.bloom \
		$(tmpdir)/e$e_l2_window1.bloom $(tmpdir)/e$e_l2_window2.bloom
	cmp $(tmpdir)/e$e_l2.bloom $(tmpdir)/e$e_l2_concat.bloom
	@echo '------------------'
	@echo '$@: PASSED'
	@echo '------------------'

#------------------------------------------------------------
# parallel_load_2_files_test
#------------------------------------------------------------

parallel_load_2_files_test: $(tmpdir)/e$e_l2.bloom \
				$(tmpdir)/e$e_1.fq \
				$(tmpdir)/e$e_2.fq
	echo 'b_div_2: $(b_div_2)'
	$(bloom) build -v -k$k -b$(b_div_2) $(tmpdir)/e$e_l1_read1.bloom \
		$(tmpdir)/e$e_1.fq
	$(bloom) build -v -k$k -b$(b_div_2) $(tmpdir)/e$e_l1_read2.bloom \
		$(tmpdir)/e$e_2.fq
	$(bloom) build -v -k$k -b$b -l2 -L 1=$(tmpdir)/e$e_l1_read2.bloom \
		$(tmpdir)/e$e_l2_read1.bloom $(tmpdir)/e$e_1.fq
	$(bloom) build -v -k$k -b$b -l2 -L 1=$(tmpdir)/e$e_l1_read1.bloom \
		$(tmpdir)/e$e_l2_read2.bloom $(tmpdir)/e$e_2.fq
	$(bloom) union -v -k$k $(tmpdir)/e$e_l2_union.bloom \
		$(tmpdir)/e$e_l2_read1.bloom \
		$(tmpdir)/e$e_l2_read2.bloom
	cmp $(tmpdir)/e$e_l2.bloom $(tmpdir)/e$e_l2_union.bloom
	@echo '------------------'
	@echo '$@: PASSED'
	@echo '------------------'

#------------------------------------------------------------
# parallel_load_3_files_test
#------------------------------------------------------------

parallel_load_3_files_test: $(tmpdir)/e$e_l2.bloom \
		$(tmpdir)/e$e_reads_1of3.fq \
		$(tmpdir)/e$e_reads_2of3.fq \
		$(tmpdir)/e$e_reads_3of3.fq
	$(bloom) build -v -k$k -b$(b_div_2) $(tmpdir)/e$e_l1_1of3.bloom \
		$(tmpdir)/e$e_reads_1of3.fq
	$(bloom) build -v -k$k -b$(b_div_2) $(tmpdir)/e$e_l1_2of3.bloom \
		$(tmpdir)/e$e_reads_2of3.fq
	$(bloom) build -v -k$k -b$(b_div_2) $(tmpdir)/e$e_l1_3of3.bloom \
		$(tmpdir)/e$e_reads_3of3.fq
	$(bloom) build -v -k$k -b$b -l2 \
		-L 1=$(tmpdir)/e$e_l1_2of3.bloom \
		-L 1=$(tmpdir)/e$e_l1_3of3.bloom \
		$(tmpdir)/e$e_l2_1of3.bloom $(tmpdir)/e$e_reads_1of3.fq
	$(bloom) build -v -k$k -b$b -l2 \
		-L 1=$(tmpdir)/e$e_l1_1of3.bloom \
		-L 1=$(tmpdir)/e$e_l1_3of3.bloom \
		$(tmpdir)/e$e_l2_2of3.bloom $(tmpdir)/e$e_reads_2of3.fq
	$(bloom) build -v -k$k -b$b -l2 \
		-L 1=$(tmpdir)/e$e_l1_1of3.bloom \
		-L 1=$(tmpdir)/e$e_l1_2of3.bloom \
		$(tmpdir)/e$e_l2_3of3.bloom $(tmpdir)/e$e_reads_3of3.fq
	$(bloom) union -v -k$k $(tmpdir)/e$e_l2_union.bloom \
		$(tmpdir)/e$e_l2_1of3.bloom \
		$(tmpdir)/e$e_l2_2of3.bloom \
		$(tmpdir)/e$e_l2_3of3.bloom
	cmp $(tmpdir)/e$e_l2.bloom $(tmpdir)/e$e_l2_union.bloom
	@echo '------------------'
	@echo '$@: PASSED'
	@echo '------------------'

#------------------------------------------------------------
# abyss_bloom_dist_1_file_test
#------------------------------------------------------------

abyss_bloom_dist_1_file_test: $(tmpdir)/e$e_l2.bloom \
		$(tmpdir)/e$e_1.fq \
		$(tmpdir)/e$e_2.fq
	cat $(tmpdir)/e$e_1.fq $(tmpdir)/e$e_2.fq > $(tmpdir)/e$e_1cat2.fq
	$(bloom_dist) -C $(tmpdir) name=dist-1-file k=$k b=$(b_div_2) w=2 \
		files='e$e_1cat2.fq'
	cmp $(tmpdir)/e$e_l2.bloom <(gunzip -c $(tmpdir)/dist-1-file.bloom.gz)
	@echo '------------------'
	@echo '$@: PASSED'
	@echo '------------------'

#------------------------------------------------------------
# abyss_bloom_dist_2_files_test
#------------------------------------------------------------

abyss_bloom_dist_2_files_test: $(tmpdir)/e$e_l2.bloom \
		$(tmpdir)/e$e_1.fq \
		$(tmpdir)/e$e_2.fq
	$(bloom_dist) -C $(tmpdir) name=dist-2-files k=$k b=$(b_div_2) w=2 \
		files='e$e_1.fq e$e_2.fq'
	cmp $(tmpdir)/e$e_l2.bloom <(gunzip -c $(tmpdir)/dist-2-files.bloom.gz)
	@echo '------------------'
	@echo '$@: PASSED'
	@echo '------------------'

#------------------------------------------------------------
# abyss_bloom_dist_3_files_test
#------------------------------------------------------------

abyss_bloom_dist_3_files_test: $(tmpdir)/e$e_l2.bloom \
		$(tmpdir)/e$e_reads_1of3.fq \
		$(tmpdir)/e$e_reads_2of3.fq \
		$(tmpdir)/e$e_reads_3of3.fq
	$(bloom_dist) -C $(tmpdir) name=dist-3-files k=$k b=$(b_div_2) w=2 \
		files='e$e_reads_1of3.fq e$e_reads_2of3.fq e$e_reads_3of3.fq'
	cmp $(tmpdir)/e$e_l2.bloom <(gunzip -c $(tmpdir)/dist-3-files.bloom.gz)
	@echo '------------------'
	@echo '$@: PASSED'
	@echo '------------------'

#------------------------------------------------------------
# abyss_bloom_illegal_chars_test
#------------------------------------------------------------

abyss_bloom_illegal_chars_test:
	$(bloom) build -v -k3 -b1M $(tmpdir)/illegal_char_test.bloom <(echo -e ">test\nAGCTagctAGCTnqrsAGCTNQRS") 
	@echo '------------------'
	@echo '$@: PASSED'
	@echo '------------------'
