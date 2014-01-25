#!/usr/bin/make -rRf

#------------------------------------------------------------
# optional params
#------------------------------------------------------------

# Number of levels in counting bloom filter
# Note: Changing this value will break the script. (I could 
# not figure out a way to write generalize the rules to work
# for any value of 'l').
l=2

#------------------------------------------------------------
# global vars
#------------------------------------------------------------

b_times_l:=$(shell echo '$b * $l' | bc)
b_mod_l:=$(shell echo '$b % $l' | bc)
l1_w1_bloom_files:=$(patsubst %,$(name)_l1_w1_%.bloom,$(files))
l1_w2_bloom_files:=$(patsubst %,$(name)_l1_w2_%.bloom,$(files))
l2_w1_bloom_files:=$(patsubst %,$(name)_l2_w1_%.bloom,$(files))
l2_w2_bloom_files:=$(patsubst %,$(name)_l2_w2_%.bloom,$(files))

PROGRAM_NAME=abyss-bloom-dist.mk

#------------------------------------------------------------
# special rules
#------------------------------------------------------------

.PHONY: args_check build
default: build

#------------------------------------------------------------
# parameter checking rule
#------------------------------------------------------------

args_check:
ifndef name
	$(error $(PROGRAM_NAME): missing required parameter 'name')
endif
ifndef k
	$(error $(PROGRAM_NAME): missing required parameter 'k')
endif
ifndef b
	$(error $(PROGRAM_NAME): missing required parameter 'b')
endif
ifndef files 
	$(error $(PROGRAM_NAME): missing required parameter 'files')
endif
ifneq ($(b_mod_l), 0)
	$(error $(PROGRAM_NAME): `b' ($b) must be divisible by `l' ($l))
endif

#------------------------------------------------------------
# main rules
#------------------------------------------------------------

build: args_check $(name).bloom

# % is one of the filenames from 'files'
$(name)_l1_w1_%.bloom: %
	abyss-bloom build -v -k$k -b$b -w1/2 $@ $<

# % is one of the filenames from 'files'
$(name)_l1_w2_%.bloom: %
	abyss-bloom build -v -k$k -b$b -w2/2 $@ $<

# % is one of the filenames from 'files'
$(name)_l2_w1_%.bloom: % $(l1_w1_bloom_files)
	abyss-bloom build -v -k$k -b$(b_times_l) -l2 -w1/2 \
		$(patsubst %,-L1=%, $(filter-out $(name)_l1_w1_$*.bloom,$(l1_w1_bloom_files))) \
		$@ $*

# % is one of the filenames from 'files'
$(name)_l2_w2_%.bloom: % $(l1_w2_bloom_files)
	abyss-bloom build -v -k$k -b$(b_times_l) -l2 -w2/2 \
		$(patsubst %,-L1=%, $(filter-out $(name)_l1_w2_$*.bloom,$(l1_w2_bloom_files))) \
		$@ $*

# final output file
$(name).bloom: $(l2_w1_bloom_files) $(l2_w2_bloom_files)
	abyss-bloom union -v -k$k $@ $^

#------------------------------------------------------------
# debugging rules
#------------------------------------------------------------

debug:
	@echo 'b_times_l=$(b_times_l)'
	@echo 'b_mod_l=$(b_mod_l)'
	@echo 'l1_w1_bloom_files="$(l1_w1_bloom_files)"'
	@echo 'l1_w2_bloom_files="$(l1_w2_bloom_files)"'
	@echo 'l2_w1_bloom_files="$(l2_w1_bloom_files)"'
	@echo 'l2_w2_bloom_files="$(l2_w2_bloom_files)"'
