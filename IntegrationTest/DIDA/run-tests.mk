#!/usr/bin/make -rRf

SHELL=/bin/bash -o pipefail

#------------------------------------------------------------
# test input/output files
#------------------------------------------------------------

# target seq for alignments
ref_url:=http://gage.cbcb.umd.edu/data/Staphylococcus_aureus/Data.original/genome.fasta
ref:=ref.fa
test_ref=test_ref.fa

# query seqs for alignments
reads_url:=http://gage.cbcb.umd.edu/data/Staphylococcus_aureus/Data.original/frag_1.fastq.gz
reads=reads.fq.gz
test_reads=test_reads.fq

# output alignment files
dida_wrapper_sam=dida_wrapper.sam
abyss_map_sam=abyss_map.sam

#------------------------------------------------------------
# params
#------------------------------------------------------------

# number of MPI tasks
np?=3
# number of threads per task
j?=1
# min align length
l?=20
# num of reads to align
n?=10000

#------------------------------------------------------------
# special targets
#------------------------------------------------------------

.PHONY: clean dida_wrapper_test

default: dida_wrapper_test

clean:
	rm -f $(dida_wrapper_sam) $(abyss_map_sam) ref-* *.lines

#------------------------------------------------------------
# downloading/building test input data
#------------------------------------------------------------

# download ref
$(ref):
	curl $(ref_url) > $@

# split ref into chunks of 100,000bp or less
$(test_ref): $(ref)
	fold -w 100000 $^ | awk '{print ">"i++; print $$0}' > $@

# download some reads
$(reads):
	curl $(reads_url) > $@

# extract first $n reads
$(test_reads): $(reads)
	zcat $(reads) | paste - - - - | head -$n | \
		tr '\t' '\n' > $@

#------------------------------------------------------------
# running DIDA/abyss-map
#------------------------------------------------------------

$(dida_wrapper_sam): $(test_reads) $(test_ref)
	abyss-dida-wrapper $^ > $@

$(abyss_map_sam): $(test_reads) $(test_ref)
	abyss-map --order -l$l $^ > $@

#------------------------------------------------------------
# tests
#------------------------------------------------------------

dida_wrapper_test: $(abyss_map_sam) $(dida_wrapper_sam)
	compare-sam $(abyss_map_sam) $(dida_wrapper_sam)
	@echo $@": PASSED!"
