---
title: konnector
author: Ben Vandervalk, Shaun Jackman, Tony Raymond, Hamid Mohamadi, Justin Chu
date: 2015-06-30
header: ABySS
footer: ABySS
section: 1
---

NAME
====

konnector - merged paired-end sequences by finding connecting paths in the de Bruijn graph

SYNOPSIS
========

`konnector -k <kmer_size> -o <output_prefix> [options]... <FASTQ> [FASTQ]...`

DESCRIPTION
===========

Konnector generates long pseudo-reads by finding connecting paths between paired-end reads within the de Bruijn graph. This can be thought of as a targeted de novo assembly in the neighbourhood of the paired-end reads. An additional feature of Konnector is the ability to extend the pseudo-reads and unmerged reads outwards, until it encounters a branching point or dead end in the de Bruijn graph.

Konnector uses a Bloom filter representation of the de Bruijn graph to minimize memory requirements, as described in: Chikhi, Rayan, and Guillaume Rizk. "Space-efficient and exact de Bruijn graph representation based on a Bloom filter." Algorithms for Molecular Biology 8.22 (2013):1.

OPTIONS
=======

Required Options
----------------
```
-k, --kmer=N               the size of a k-mer [required]
-o, --output-prefix=FILE   prefix of output FASTA files [required]
```

Bloom Filter Options
--------------------
```
-b, --bloom-size=N         size of bloom filter [500M]
-c, --min-coverage=N       kmer coverage threshold for error correction [2].
                           This option specifies the number of levels in the
                           cascading Bloom filter; it has no effect if the Bloom
                           filter is loaded from an external file.
-i, --input-bloom=FILE     load bloom filter from FILE; Bloom filter files can
                           be created separately with the 'abyss-bloom' program
```

Graph Search Limits
-------------------
```
-B, --max-branches=N       max branches in de Bruijn graph traversal;
                           use 'nolimit' for no limit [350]
-f, --min-frag=N           min fragment size in base pairs [0]
-F, --max-frag=N           max fragment size in base pairs [1000]
-P, --max-paths=N          merge at most N alternate paths; use 'nolimit'
                           for no limit [2]
```

Sequence Identity Limits
------------------------
```
-m, --read-mismatches=N    max mismatches between paths and reads; use
                           'nolimit' for no limit [nolimit]
-M, --max-mismatches=N     max mismatches between all alternate paths;
                           use 'nolimit' for no limit [2]
-x, --read-identity=N      min percent seq identity between consensus seq
                           and reads [0]
-X, --path-identity=N      min percent seq identity across alternate
                           connecting paths [0]
```

Input Options
-------------
```
-q, --trim-quality=N       trim bases from the ends of reads whose
                           quality is less than the threshold
    --standard-quality     zero quality is `!' (33), typically
                           for FASTQ and SAM files [default]
    --illumina-quality     zero quality is `@' (64), typically
                           for qseq and export files
    --chastity             discard unchaste reads [default]
    --no-chastity          do not discard unchaste reads
    --trim-masked          trim masked bases from the ends of reads
    --no-trim-masked       do not trim masked bases from the ends
                           of reads [default]
-I, --interleaved          input reads files are interleaved
```

Output Options
--------------
```
    --fastq                output merged reads in FASTQ format
                           (default is FASTA); bases that are corrected
						   or inserted by konnector are assigned a
						   fixed quality score determined by -Q
-Q, --corrected-qual       quality score for bases corrected or inserted
                           by konnector; only relevant when --fastq is
                           in effect [40]
    --mask                 mask new and changed bases as lower case
    --no-mask              do not mask bases [default]
-p, --alt-paths-mode       output a separate pseudoread for each alternate
                           path connecting a read pair (default is to create
                           a consensus sequence of all connecting paths).
						   The limit on the number of alternate paths is
						   specified by the '--max-paths' option.
                           The sequence IDs for alternate paths are named:
						   ${orig_read_id}_1, ${orig_read_id}_2, ...
--preserve-reads           don't correct any bases within the reads [disabled]
-v, --verbose              display verbose output
```

Debugging Options
-----------------
```
-d, --dot-file=FILE        write graph traversals to a DOT file
-r, --read-name=STR        only process reads with names that contain STR
-t, --trace-file=FILE      write graph search stats to FILE
```

Sequence Extension Options
--------------------------
```
-D, --dup-bloom-size=N     use an additional Bloom filter to avoid
                           assembling the same region of the genome
                           multiple times. This option is highly
                           recommended whenever -E (--extend) is used
                           and has no effect otherwise. As a rule of
                           thumb, the Bloom filter size should be
                           about twice the target genome size [disabled]
-E, --extend               in addition to connecting read pairs,
                           extend the merged reads outwards to the next
                           dead end or branching point in the de Brujin
                           graph. For read pairs that were not successfully
                           connected, trim the single-end reads at both ends
						   and extend them independently.
```

Other Options
-------------
```
-e, --fix-errors           find and fix single-base errors when reads
                           have no kmers in bloom filter [disabled]
-j, --threads=N            use N parallel threads [1]
-n  --no-limits            disable all limits; equivalent to
                           '-B nolimit -m nolimit -M nolimit -P nolimit'
    --help                 display this help and exit
    --version              output version information and exit
```

OUTPUT FILES
============

`$PREFIX` in the filenames below is determined by the `-o` option.

Without `--extend`:

  * `$PREFIX_pseudoreads.fa`: Pseudo-reads created by connecting paired-end reads.
  * `$PREFIX_reads_1.fq`: Read 1 from read pairs that could not be connected.
  * `$PREFIX_reads_2.fq`: Read 2 from read pairs that could not be connected.

With `--extend`:

  * `$prefix_pseudoreads.fa`: Pseudo-reads created by connecting paired-end reads, which may or may not be extended. Also contains single-end reads from read pairs that could not be connected, but which could be trimmed and/or extended.
  * `$PREFIX_reads_1.fq`: Read 1 from read pairs that could not be connected and which could not be trimmed (because they contain no "good" k-mers).
  * `$PREFIX_reads_2.fq`: Read 2 from read pairs that could not be connected and which could not be trimmed (because they contain no "good" k-mers).

AUTHORS
=======

Ben Vandervalk, Shaun Jackman, Tony Raymond, Hamid Mohamadi, Justin Chu.

REPORTING BUGS
==============

Report bugs to <abyss-users@bcgsc.ca>.
