Description
===========

Sealer is an application of Konnector that closes intra-scaffold gaps. It performs three sequential functions. First, regions with Ns are identified from an input scaffold. Flanking nucleotues (2 x 100bp) are extracted from those regions while respecting the strand (5' to 3') direction on the sequence immediately downstream of each gap. In the second step, flanking sequence pairs are used as input to Konnector along with a set of reads with a high level of coverage redundancy. Ideally, the reads should represent the original dataset from which the draft assembly is generated, or further whole genome shotgun (WGS) sequencing data generated from the same sample. Within Konnector, the input WGS reads are used to populate a Bloom filter, tiling the reads with a sliding window of length k, thus generating a probabilistic representation of all the k-mers in the reads. Konnector also uses crude error removal and correctional algorithms, eliminating singletons (k-mers that are observed only once) and fixing base mismatches in the flanking sequence pairs. Sealer launches Konnector processes using a user-input range of k-mer lengths. In the third and final operation, succesfully merged sequences are inserted into the gaps of the original scaffolds, and Sealer outputs a new gap-filled scaffold file. 

Installation
============

See ABySS installation instructions.


How to run as stand-alone application
=====================================

`abyss-sealer [-k values...] [-o outputprefix] [-S assembly file] [options...] [reads...]`

Sealer requires the following information to run:
- draft assembly
- user-supplied k values (>0)
- output prefix
- WGS reads  (for building Bloom Filters) 


Sample commands
===============

Without pre-built bloom filters:

`abyss-sealer -k90 -k80 -o run1 -S test.fa read1.fa.gz read2.fa.gz`

With pre-built bloom filters:

`abyss-sealer -k90 -k80 -o run1 -S test.fa -i k90.bloom -i k80.bloom read1.fa.gz read2.fa.gz`


Suggested parameters for first run
==================================

When running Sealer for the first time on a dataset, we recommend using the following parameters. *P* is the threshold for number of paths allowed to be traversed. When set to 10, Konnector will attempt to close gaps even when there are 10 different paths found. It would attempt to create a consensus sequence between these paths. The default setting is 2. The required parameter *k* sets the k-mer length, and these suggested values will provide insight for optimizing Sealer. For instance, if after the first run most gaps are closed at k90, then on the next run, using k92, k91, k90, k89, and k88 may close even more gaps. Note that these suggested k values are only applicable for datasets without a pre-built Bloom Filter. Note: the max possible k value for Sealer is 100. 

- `-P 10` 
- `-k90 -k80 -k70 -k60 -k50 -k40 -k30`


Output files
============

- prefix_log.txt
- prefix_scaffold.fa
- prefix_merged.fa
- prefix_flanks_1.fq  -> if `--print-flanks` option used
- prefix_flanks_2.fq  -> if `--print-flanks` option used

The log file contains results of each Konnector run. The structure of one run is as follows:

* \#\# unique gaps closed for k## 	< # closed gaps with unique path + # closed gaps with multiple paths
* No start/goal kmer: ###			< # unclosed gaps with no start/goal k-mer
* No path: ###				< # unclosed gaps with no path found
* Unique path: ###			< # gaps that closed with unique paths
* Multiple paths: ###			< # gaps that closed with >1 paths. A consensus sequence between the paths are output. This will likely contain ambiguity codes
* Too many paths: ###			< # unclosed gaps with paths greater than allowed by -P parameter (default 2)
* Too many branches: ###			< # unclosed gaps with branches greater than allowed by -B paramter (default 1000)
* Too many path/path mismatches: ###	< # unclosed gaps with path/path mismatches greater than allowed by -M paramter (default no-limit)
* Too many path/read mismatches: ###	< # unclosed gaps with path/read mismatches greater than allowed by -m parameter (default no-limit)
* Contains cycle: ###			< # unclosed gaps containing cycles
* Exceeded mem limit: ###			
* Skipped: ###				< # gaps skipped
* \#\#\# flanks left				< # gaps left unclosed and will be inserted to next K run. Closed gaps no longer inserted into subsequent K runs.
* k## run complete
* Total gaps closed so far = ###		< # gaps closed by all K runs so far

The scaffold.fa file is a gap-filled version of the draft assembly inserted into Sealer. The merged.fa file contains every newly generated sequence that were inserted into gaps, including the flanking sequences. Negative sizes of new sequences indicate Konnector collapsed the pair of flanking sequences. For example:

\>[scaffold ID]\_[original start position of gap on scaffold]\_[size of new sequence]
ACGCGACGAGCAGCGAGCACGAGCAGCGACGAGCGACGACGAGCAGCGACGAGCG


If `--print-flanks` option is enabled, Sealer outputs the flanking sequences
used to insert into Konnector. This may be useful should users which to double
check if this tool is extracting the correct sequences surrounding gaps. The
structure of these files are as follows:

\>[scaffold ID]\_[original start position of gap on scaffold]\_[size of gap]/[1 or 2 indicating whether left or right flank]
GCTAGCTAGCTAGCTGATCGATCGTAGCTAGCTAGCTGACTAGCTGATCAGTCGA


How to optimize for gap closure
===============================

To optimize Sealer, users can observe the log files generated after a run and
adjust parameters accordingly. If K runs are showing gaps having too many
paths or branches, consider increasing -P or -B parameters, respectively. 

Also consider increasing the number of K values used. Generally, large K-mers
are better able to address highly repetitive genomic regions, while smaller
K-mers are better able to resolve areas of low coverage. 

Sometimes, datasets will only have gaps closed around a certain k value.
After running Sealer with the suggested k parameters of -k90 to -k30 (interval
of 10), observe which k value has the most gap closures. Supposing gaps are
mostly being closed at k90, then consider running Sealer with k values around
k90. i.e. -k95 to k85 (interval of 1)


Runtime and memory usage
========================

More K values mean more bloom filters will be required, which will increase
runtime as it takes time to build/load each bloom filter at the beginning of
each k run. Memory usage is not affected by using more bloom filters. 

The larger value used for parameters such as `-P`, `-B` or `-F` will increase
runtime. 


Options
=======

Parameters of `abyss-sealer`

* `--print-flanks`: outputs flank files
* `-S`,`--input-scaffold=FILE`:	load scaffold from FILE
* `-L`,`--flank-length=N`: length of flanks to be used as pseudoreads [`100`]
* `-D`,`--flank-distance=N`: distance of flank from gap [0]
* `-j`,`--threads=N`: use N parallel threads [1]
* `-k`,`--kmer=N`: the size of a k-mer
* `-b`,`--bloom-size=N`: size of bloom filter [500M]
* `-B`,`--max-branches=N`: max branches in de Bruijn graph traversal; use 'nolimit' for no limit [1000]
* `-d`,`--dot-file=FILE`: write graph traversals to a DOT file
* `-e`,`--fix-errors`: find and fix single-base errors when reads have no kmers in bloom filter [disabled]
* `-f`,`--min-frag=N`: min fragment size in base pairs [0]
* `-F`,`--max-frag=N`: max fragment size in base pairs [1000]
* `-i`,`--input-bloom=FILE`: load bloom filter from FILE
* `--mask`: mask new and changed bases as lower case
* `--no-mask`: do not mask bases [default]
* `--chastity`: discard unchaste reads [default]
* `--no-chastity`: do not discard unchaste reads
* `--trim-masked`: trim masked bases from the ends of reads
* `--no-trim-masked`: do not trim masked bases from the ends of reads [default]
* `-l`,`--long-search`: start path search as close as possible to the beginnings of reads. Takes more time but improves results when bloom filter false positive rate is high [disabled]
* `-m,`--flank-mismatches=N`: max mismatches between paths and flanks; use 'nolimit' for no limit [nolimit]
* `-M,`--max-mismatches=N`: max mismatches between all alternate paths; use 'nolimit' for no limit [nolimit]
* `-n `--no-limits`: disable all limits; equivalent to '-B nolimit -m nolimit -M nolimit -P nolimit'
* `-o,`--output-prefix=FILE`: prefix of output FASTA files [required]
* `-P,`--max-paths=N`: merge at most N alternate paths; use 'nolimit' for no limit [2]
* `-q,`--trim-quality=N`: trim bases from the ends of reads whose quality is less than the threshold
* `--standard-quality`: zero quality is `!' (33) default for FASTQ and SAM files
* `--illumina-quality`: zero quality is `@' (64) default for qseq and export files
* `-r,`--read-name=STR`: only process reads with names that contain STR
* `-s,`--search-mem=N`: mem limit for graph searches; multiply by the number of threads (-j) to get the total mem used for graph traversal [500M]
* `-t,`--trace-file=FILE`: write graph search stats to FILE
* `-v,`--verbose`: display verbose output
* `--help`: display this help and exit
* `--version`: output version information and exit


