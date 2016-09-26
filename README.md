ABySS
=====

ABySS is a *de novo* sequence assembler intended for short paired-end
reads and large genomes.

Contents
========

* [Quick Start](#quick-start)
	* [Install ABySS on Debian or Ubuntu](#install-abyss-on-debian-or-ubuntu)
	* [Install ABySS on Mac OS X](#install-abyss-on-mac-os-x)
* [Dependencies](#dependencies)
* [Compiling ABySS from GiHub](#compiling-abyss-from-github)
* [Compiling ABySS from source](#compiling-abyss-from-source)
* [Assembling a paired-end library](#assembling-a-paired-end-library)
* [Assembling multiple libraries](#assembling-multiple-libraries)
* [Scaffolding](#scaffolding)
* [Rescaffolding with long sequences](#rescaffolding-with-long-sequences)
* [Assembling using a Bloom filter de Bruijn graph](#assembling-using-a-bloom-filter-de-bruijn-graph)
* [Assembling using a paired de Bruijn graph](#assembling-using-a-paired-de-bruijn-graph)
* [Assembling a strand-specific RNA-Seq library](#assembling-a-strand-specific-rna-seq-library)
* [Optimizing the parameter k](#optimizing-the-parameter-k)
* [Parallel processing](#parallel-processing)
* [Running ABySS on a cluster](#running-abyss-on-a-cluster)
* [Using the DIDA alignment framework](#using-the-dida-alignment-framework)
* [Assembly Parameters](#assembly-parameters)
* [ABySS programs](#abyss-programs)
* [Export to SQLite Database](#export-to-sqlite-database)
* [Publications](#publications)
* [Support](#support)
* [Authors](#authors)

Quick Start
===========

## Install ABySS on Debian or Ubuntu

Run the command

	sudo apt-get install abyss

or download and install the
[Debian package](http://www.bcgsc.ca/platform/bioinfo/software/abyss).

## Install ABySS on Mac OS X

Install [Homebrew](http://brew.sh/), and run the commands

	brew install homebrew/science/abyss

## Assemble a small synthetic data set

	wget http://www.bcgsc.ca/platform/bioinfo/software/abyss/releases/1.3.4/test-data.tar.gz
	tar xzvf test-data.tar.gz
	abyss-pe k=25 name=test \
		in='test-data/reads1.fastq test-data/reads2.fastq'

## Calculate assembly contiguity statistics

	abyss-fac test-unitigs.fa

Dependencies
============

ABySS requires the following libraries:

* [Boost](http://www.boost.org/)
* [Open MPI](http://www.open-mpi.org)
* [sparsehash](https://code.google.com/p/sparsehash/)
* [SQLite](http://www.sqlite.org/)

ABySS requires a C++ compiler that supports
[OpenMP](http://www.openmp.org) such as [GCC](http://gcc.gnu.org).

ABySS will receive an error when compiling with Boost 1.51.0 or 1.52.0
since they contain a bug. Later versions of Boost compile without error.

Compiling ABySS from GitHub
===========================

When installing ABySS from GitHub source the following tools are
required:

* [Autoconf](http://www.gnu.org/software/autoconf)
* [Automake](http://www.gnu.org/software/automake)

To generate the configure script and make files:

	./autogen.sh

See "Compiling ABySS from source" for further steps.

Compiling ABySS from source
===========================

To compile and install ABySS in `/usr/local`:

	./configure
	make
	sudo make install

To install ABySS in a specified directory:

	./configure --prefix=/opt/abyss
	make
	sudo make install

ABySS uses OpenMP for parallelization, which requires a modern
compiler such as GCC 4.2 or greater. If you have an older compiler, it
is best to upgrade your compiler if possible. If you have multiple
versions of GCC installed, you can specify a different compiler:

	./configure CC=gcc-4.6 CXX=g++-4.6

ABySS requires the Boost C++ libraries. Many systems come with Boost
installed. If yours does not, you can download
[Boost](http://www.boost.org/users/download).
It is not necessary to compile Boost before installing it. The Boost
header file directory should be found at `/usr/include/boost`, in the
ABySS source directory, or its location specified to `configure`:

	./configure --with-boost=/usr/local/include

If you wish to build the parallel assembler with MPI support,
MPI should be found in `/usr/include` and `/usr/lib` or its location
specified to `configure`:

	./configure --with-mpi=/usr/lib/openmpi

ABySS should be built using the sparsehash library to reduce memory
usage, although it will build without. sparsehash should be found in
`/usr/include` or its location specified to `configure`:

	./configure CPPFLAGS=-I/usr/local/include

If SQLite is installed in non-default directories, its location can be
specified to `configure`:

	./configure --with-sqlite=/opt/sqlite3

The default maximum k-mer size is 64 and may be decreased to reduce
memory usage or increased at compile time. This value must be a
multiple of 32 (i.e. 32, 64, 96, 128, etc):

	./configure --enable-maxk=96

If you encounter compiler warnings, you may ignore them like so:

	make AM_CXXFLAGS=-Wall

To run ABySS, its executables should be found in your `PATH`. If you
installed ABySS in `/opt/abyss`, add `/opt/abyss/bin` to your `PATH`:

	PATH=/opt/abyss/bin:$PATH

Assembling a paired-end library
===============================

To assemble paired reads in two files named `reads1.fa` and
`reads2.fa` into contigs in a file named `ecoli-contigs.fa`, run the
command:

	abyss-pe name=ecoli k=64 in='reads1.fa reads2.fa'

The parameter `in` specifies the input files to read, which may be in
FASTA, FASTQ, qseq, export, SRA, SAM or BAM format and compressed with
gz, bz2 or xz and may be tarred. The assembled contigs will be stored
in `${name}-contigs.fa`.

A pair of reads must be named with the suffixes `/1` and `/2` to
identify the first and second read, or the reads may be named
identically. The paired reads may be in separate files or interleaved
in a single file.

Reads without mates should be placed in a file specified by the
parameter `se` (single-end). Reads without mates in the paired-end
files will slow down the paired-end assembler considerably during the
`abyss-fixmate` stage.

Assembling multiple libraries
=============================

The distribution of fragment sizes of each library is calculated
empirically by aligning paired reads to the contigs produced by the
single-end assembler, and the distribution is stored in a file with
the extension `.hist`, such as `ecoli-3.hist`. The N50 of the
single-end assembly must be well over the fragment-size to obtain an
accurate empirical distribution.

Here's an example scenario of assembling a data set with two different
fragment libraries and single-end reads. Note that the names of the libraries
(`pea` and `peb`) are arbitrary.

 * Library `pea` has reads in two files,
   `pea_1.fa` and `pea_2.fa`.
 * Library `peb` has reads in two files,
   `peb_1.fa` and `peb_2.fa`.
 * Single-end reads are stored in two files, `se1.fa` and `se2.fa`.

The command line to assemble this example data set is:

	abyss-pe k=64 name=ecoli lib='pea peb' \
		pea='pea_1.fa pea_2.fa' peb='peb_1.fa peb_2.fa' \
		se='se1.fa se2.fa'

The empirical distribution of fragment sizes will be stored in two
files named `pea-3.hist` and `peb-3.hist`. These files may be
plotted to check that the empirical distribution agrees with the
expected distribution. The assembled contigs will be stored in
`${name}-contigs.fa`.

Scaffolding
===========

Long-distance mate-pair libraries may be used to scaffold an assembly.
Specify the names of the mate-pair libraries using the parameter `mp`.
The scaffolds will be stored in the file `${name}-scaffolds.fa`.
Here's an example of assembling a data set with two paired-end
libraries and two mate-pair libraries. Note that the names of the libraries
(`pea`, `peb`, `mpa`, `mpb`) are arbitrary.

	abyss-pe k=64 name=ecoli lib='pea peb' mp='mpc mpd' \
		pea='pea_1.fa pea_2.fa' peb='peb_1.fa peb_2.fa' \
		mpc='mpc_1.fa mpc_2.fa' mpd='mpd_1.fa mpd_2.fa'

The mate-pair libraries are used only for scaffolding and do not
contribute towards the consensus sequence.

Rescaffolding with long sequences
=================================

Long sequences such as RNA-Seq contigs can be used to rescaffold an
assembly. Sequences are aligned using BWA-MEM to the assembled
scaffolds. Additional scaffolds are then formed between scaffolds that
can be linked unambiguously when considering all BWA-MEM alignments.

Similar to scaffolding, the names of the datasets can be specified with
the `long` parameter. These scaffolds will be stored in the file
`${name}-trans-scaffs.fa`. The following is an example of an assembly with PET, MPET and an RNA-Seq assembly. Note that the names of the libraries are arbitrary.

	abyss-pe k=64 name=ecoli lib='pe1 pe2' mp='mp1 mp2' long='longa' \
		pe1='pe1_1.fa pe1_2.fa' pe2='pe2_1.fa pe2_2.fa' \
		mp1='mp1_1.fa mp1_2.fa' mp2='mp2_1.fa mp2_2.fa' \
		longa='longa.fa'

Assembling using a Bloom filter de Bruijn graph
=========================================

Assemblies may be performed using a _Bloom filter de Bruijn graph_, which
typically reduces memory requirements by an order of magnitude. To assemble in
Bloom filter mode, the user must specify 3 additional parameters: `B` (Bloom
filter size in bytes), `H` (number of Bloom filter hash functions), and `kc`
(minimum k-mer count threshold). `B` is the overall memory budget for the Bloom
filter assembler, and may be specified with unit suffixes 'k' (kilobytes), 'M'
(megabytes), 'G' (gigabytes). If no units are specified bytes are assumed. For
example, the following will run a E. coli assembly with an overall memory budget
of 100 megabytes, 3 hash functions, a minimum k-mer count threshold of 3, with
verbose logging enabled:

	abyss-pe name=ecoli k=64 in='reads1.fa reads2.fa' B=100M H=3 kc=3 v=-v

At the current time, the user must calculate suitable values for `B` and `H` on
their own, and finding the best value for `kc` may require experimentation
(optimal values are typically in the range of 2-4). Internally, the Bloom filter
assembler divides the memory budget (`B`) equally across (`kc` + 1) Bloom
filters, where `kc` Bloom filters are used for the cascading Bloom filter and
one additional Bloom filter is used to track k-mers that have previously been
included in contigs. Users are recommended to target a Bloom filter false
positive rate (FPR) that is less than 5%, as reported by the assembly log when
using the `v=-v` option (verbose level 1).

Assembling using a paired de Bruijn graph
=========================================

Assemblies may be performed using a _paired de Bruijn graph_ instead
of a standard de Bruijn graph.  In paired de Bruijn graph mode, ABySS
uses _k-mer pairs_ in place of k-mers, where each k-mer pair consists of
two equal-size k-mers separated by a fixed distance.  A k-mer pair
is functionally similar to a large k-mer spanning the breadth of the k-mer
pair, but uses less memory because the sequence in the gap is not stored.
To assemble using paired de Bruijn graph mode, specify both individual
k-mer size (`K`) and k-mer pair span (`k`). For example, to assemble E.
coli with a individual k-mer size of 16 and a k-mer pair span of 64:

	abyss-pe name=ecoli K=16 k=64 in='reads1.fa reads2.fa'

In this example, the size of the intervening gap between k-mer pairs is
32 bp (64 - 2\*16). Note that the `k` parameter takes on a new meaning
in paired de Bruijn graph mode. `k` indicates kmer pair span in
paired de Bruijn graph mode (when `K` is set), whereas `k` indicates
k-mer size in standard de Bruijn graph mode (when `K` is not set).

Assembling a strand-specific RNA-Seq library
============================================

Strand-specific RNA-Seq libraries can be assembled such that the
resulting unitigs, contigs and scaffolds are oriented correctly with
respect to the original transcripts that were sequenced. In order to
run ABySS in strand-specific mode, the `SS` parameter must be used as
in the following example:

	abyss-pe name=SS-RNA k=64 in='reads1.fa reads2.fa' SS=--SS

The expected orientation for the read sequences with respect to the
original RNA is RF. i.e. the first read in a read pair is always in
reverse orientation.

Optimizing the parameter k
==========================

To find the optimal value of `k`, run multiple assemblies and inspect
the assembly contiguity statistics. The following shell snippet will
assemble for every eighth value of `k` from 50 to 90.

	for k in `seq 50 8 90`; do
		mkdir k$k
		abyss-pe -C k$k name=ecoli k=$k in=../reads.fa
	done
	abyss-fac k*/ecoli-contigs.fa

The default maximum value for `k` is 96. This limit may be changed at
compile time using the `--enable-maxk` option of configure. It may be
decreased to 32 to decrease memory usage or increased to larger values.

Parallel processing
===================

The `np` option of `abyss-pe` specifies the number of processes to
use for the parallel MPI job. Without any MPI configuration, this will
allow you to use multiple cores on a single machine. To use multiple
machines for assembly, you must create a `hostfile` for `mpirun`,
which is described in the `mpirun` man page.

*Do not* run `mpirun -np 8 abyss-pe`. To run ABySS with 8 threads, use
`abyss-pe np=8`. The `abyss-pe` driver script will start the MPI
process, like so: `mpirun -np 8 ABYSS-P`.

The paired-end assembly stage is multithreaded, but must run on a
single machine. The number of threads to use may be specified with the
parameter `j`. The default value for `j` is the value of `np`.

Running ABySS on a cluster
==========================

ABySS integrates well with cluster job schedulers, such as:

 * SGE (Sun Grid Engine)
 * Portable Batch System (PBS)
 * Load Sharing Facility (LSF)
 * IBM LoadLeveler

For example, to submit an array of jobs to assemble every eighth value of
`k` between 50 and 90 using 64 processes for each job:

	qsub -N ecoli -pe openmpi 64 -t 50-90:8 \
		<<<'mkdir k$SGE_TASK_ID && abyss-pe -C k$SGE_TASK_ID in=/data/reads.fa'

Using the DIDA alignment framework
=================================

ABySS supports the use of DIDA (Distributed Indexing Dispatched Alignment),
an MPI-based framework for computing sequence alignments in parallel across
multiple machines. The DIDA software must be separately downloaded and
installed from http://www.bcgsc.ca/platform/bioinfo/software/dida. In
comparison to the standard ABySS alignment stages which are constrained
to a single machine, DIDA offers improved performance and the ability to
scale to larger targets. Please see the DIDA section of the abyss-pe man
page (in the `doc` subdirectory) for details on usage.

Assembly Parameters
===================

Parameters of the driver script, `abyss-pe`

 * `a`: maximum number of branches of a bubble [`2`]
 * `b`: maximum length of a bubble (bp) [`""`]
 * `B`: Bloom filter size (e.g. "100M")
 * `c`: minimum mean k-mer coverage of a unitig [`sqrt(median)`]
 * `d`: allowable error of a distance estimate (bp) [`6`]
 * `e`: minimum erosion k-mer coverage [`round(sqrt(median))`]
 * `E`: minimum erosion k-mer coverage per strand [1 if sqrt(median) > 2 else 0]
 * `G`: genome size, used to calculate NG50 [disabled]
 * `H`: number of Bloom filter hash functions [1]
 * `j`: number of threads [`2`]
 * `k`: size of k-mer (when `K` is not set) or the span of a k-mer pair (when `K` is set)
 * `kc`: minimum k-mer count threshold for Bloom filter assembly [`2`]
 * `K`: the length of a single k-mer in a k-mer pair (bp)
 * `l`: minimum alignment length of a read (bp) [`40`]
 * `m`: minimum overlap of two unitigs (bp) [`30`]
 * `n`: minimum number of pairs required for building contigs [`10`]
 * `N`: minimum number of pairs required for building scaffolds [`n`]
 * `np`: number of MPI processes [`1`]
 * `p`: minimum sequence identity of a bubble [`0.9`]
 * `q`: minimum base quality [`3`]
 * `s`: minimum unitig size required for building contigs (bp) [`1000`]
 * `S`: minimum contig size required for building scaffolds (bp) [`1000-10000`]
 * `t`: maximum length of blunt contigs to trim [`k`]
 * `v`: use `v=-v` for verbose logging, `v=-vv` for extra verbose [`disabled`]
 * `x`: spaced seed (Bloom filter assembly only)

Please see the
[abyss-pe](http://manpages.ubuntu.com/abyss-pe.1.html)
manual page for more information on assembly parameters.

Environment variables
=====================

`abyss-pe` configuration variables may be set on the command line or from the environment, for example with `export k=20`. It can happen that `abyss-pe` picks up such variables from your environment that you had not intended, and that can cause trouble. To troubleshoot that situation, use the `abyss-pe env` command to print the values of all the `abyss-pe` configuration variables:

	abyss-pe env [options]

ABySS programs
==============

`abyss-pe` is a driver script implemented as a Makefile. Any option of
`make` may be used with `abyss-pe`. Particularly useful options are:

 * `-C dir`, `--directory=dir`
   Change to the directory `dir` and store the results there.
 * `-n`, `--dry-run`
   Print the commands that would be executed, but do not execute
   them.

`abyss-pe` uses the following programs, which must be found in your
`PATH`:

 * `ABYSS`: de Bruijn graph assembler
 * `ABYSS-P`: parallel (MPI) de Bruijn graph assembler
 * `AdjList`: find overlapping sequences
 * `DistanceEst`: estimate the distance between sequences
 * `MergeContigs`: merge sequences
 * `MergePaths`: merge overlapping paths
 * `Overlap`: find overlapping sequences using paired-end reads
 * `PathConsensus`: find a consensus sequence of ambiguous paths
 * `PathOverlap`: find overlapping paths
 * `PopBubbles`: remove bubbles from the sequence overlap graph
 * `SimpleGraph`: find paths through the overlap graph
 * `abyss-fac`: calculate assembly contiguity statistics
 * `abyss-filtergraph`: remove shim contigs from the overlap graph
 * `abyss-fixmate`: fill the paired-end fields of SAM alignments
 * `abyss-map`: map reads to a reference sequence
 * `abyss-scaffold`: scaffold contigs using distance estimates
 * `abyss-todot`: convert graph formats and merge graphs

This [flowchart](https://github.com/bcgsc/abyss/blob/master/doc/flowchart.pdf) shows the ABySS assembly pipeline its intermediate files.

Export to SQLite Database
=========================

ABySS has a built-in support for SQLite database to export log values into a SQLite file and/or `.csv` files at runtime.

## Database parameters
Of `abyss-pe`:
 * `db`: path to SQLite repository file [`$(name).sqlite`]
 * `species`: name of species to archive [ ]
 * `strain`: name of strain to archive [ ]
 * `library`: name of library to archive [ ]

For example, to export data of species 'Ecoli', strain 'O121' and library 'pea' into your SQLite database repository named '/abyss/test.sqlite':

	abyss-pe db=/abyss/test.sqlite species=Ecoli strain=O121 library=pea [other options]

## Helper programs
Found in your `path`:
 * `abyss-db-txt`: create a flat file showing entire repository at a glance
 * `abyss-db-csv`: create `.csv` table(s) from the repository

Usage:

    abyss-db-txt /your/repository
    abyss-db-csv /your/repository program(s)

For example,

	abyss-db-txt repo.sqlite

	abyss-db-csv repo.sqlite DistanceEst
	abyss-db-csv repo.sqlite DistanceEst abyss-scaffold
	abyss-db-csv repo.sqlite --all

Publications
============

## [ABySS](http://genome.cshlp.org/content/19/6/1117)

Simpson, Jared T., Kim Wong, Shaun D. Jackman, Jacqueline E. Schein,
Steven JM Jones, and Inanc Birol.
**ABySS: a parallel assembler for short read sequence data**.
*Genome research* 19, no. 6 (2009): 1117-1123.
[doi:10.1101/gr.089532.108](http://dx.doi.org/10.1101/gr.089532.108)

## [Trans-ABySS](http://www.nature.com/nmeth/journal/v7/n11/abs/nmeth.1517.html)

Robertson, Gordon, Jacqueline Schein, Readman Chiu, Richard Corbett,
Matthew Field, Shaun D. Jackman, Karen Mungall et al.
**De novo assembly and analysis of RNA-seq data**.
*Nature methods* 7, no. 11 (2010): 909-912.
[doi:10.1038/10.1038/nmeth.1517](http://dx.doi.org/10.1038/nmeth.1517)

## [ABySS-Explorer](http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=5290690)

Nielsen, Cydney B., Shaun D. Jackman, Inanc Birol, and Steven JM Jones.
**ABySS-Explorer: visualizing genome sequence assemblies**.
*IEEE Transactions on Visualization and Computer Graphics*
15, no. 6 (2009): 881-888.
[doi:10.1109/TVCG.2009.116](http://dx.doi.org/10.1109/TVCG.2009.116)

Support
=======

[Ask a question](https://www.biostars.org/p/new/post/?tag_val=abyss,assembly)
on [Biostars](https://www.biostars.org/t/abyss/).

[Create a new issue](https://github.com/bcgsc/abyss/issues) on GitHub.

Subscribe to the
[ABySS mailing list]
(http://groups.google.com/group/abyss-users),
<abyss-users@googlegroups.com>.

For questions related to transcriptome assembly, contact the
[Trans-ABySS mailing list]
(http://groups.google.com/group/trans-abyss),
<trans-abyss@googlegroups.com>.

Authors
=======

+ **[Shaun Jackman](http://sjackman.ca)** - [GitHub/sjackman](https://github.com/sjackman) - [@sjackman](https://twitter.com/sjackman)
+ **Tony Raymond** - [GitHub/traymond](https://github.com/traymond)
+ **Ben Vandervalk** - [GitHub/benvvalk ](https://github.com/benvvalk)
+ **Jared Simpson** - [GitHub/jts](https://github.com/jts)

Supervised by [**Dr. Inanc Birol**](http://www.bcgsc.ca/faculty/inanc-birol).

Copyright 2016 Canada's Michael Smith Genome Sciences Centre
