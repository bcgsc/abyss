Title: ABySS README
Author: Shaun Jackman, Anthony Raymond
Affiliation: Canada's Michael Smith Genome Sciences Centre
CSS: README.css

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
* [Assembling a strand-specific RNA-Seq library](#assembling-a-strand-specific-rna-seq-library)
* [Optimizing the parameter k](#optimizing-the-parameter-k)
* [Parallel processing](#parallel-processing)
* [Running ABySS on a cluster](#running-abyss-on-a-cluster)
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

* [Boost](http://www.boost.org)
* [sparsehash](http://code.google.com/p/sparsehash)
* [Open MPI](http://www.open-mpi.org)

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
fragment libraries and single-end reads:

 * Library `pe200` has reads in two files,
   `pe200_1.fa` and `pe200_2.fa`.
 * Library `pe500` has reads in two files,
   `pe500_1.fa` and `pe500_2.fa`.
 * Single-end reads are stored in two files, `se1.fa` and `se2.fa`.

The command line to assemble this example data set is:

	abyss-pe k=64 name=ecoli lib='pe200 pe500' \
		pe200='pe200_1.fa pe200_2.fa' pe500='pe500_1.fa pe500_2.fa' \
		se='se1.fa se2.fa'

The empirical distribution of fragment sizes will be stored in two
files named `pe200-3.hist` and `pe500-3.hist`. These files may be
plotted to check that the empirical distribution agrees with the
expected distribution. The assembled contigs will be stored in
`${name}-contigs.fa`.

Scaffolding
===========

Long-distance mate-pair libraries may be used to scaffold an assembly.
Specify the names of the mate-pair libraries using the parameter `mp`.
The scaffolds will be stored in the file `${name}-scaffolds.fa`.
Here's an example of assembling a data set with two paired-end
libraries and two mate-pair libraries:

	abyss-pe k=64 name=ecoli lib='pe1 pe2' mp='mp1 mp2' \
		pe1='pe1_1.fa pe1_2.fa' pe2='pe2_1.fa pe2_2.fa' \
		mp1='mp1_1.fa mp1_2.fa' mp2='mp2_1.fa mp2_2.fa'

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
`${name}-trans-scaffs.fa`. The following is an example of an assembly with PET, MPET and an RNA-Seq assembly:

	abyss-pe k=64 name=ecoli lib='pe1 pe2' mp='mp1 mp2' long=long1 \
		pe1='pe1_1.fa pe1_2.fa' pe2='pe2_1.fa pe2_2.fa' \
		mp1='mp1_1.fa mp1_2.fa' mp2='mp2_1.fa mp2_2.fa' \
		long1=long1.fa

Assembling a strand-specific RNA-Seq library
============================================

Strand-specific RNA-Seq libraries can be assembled such that the
resulting unitigs, conitgs and scaffolds are oriented correctly with
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
assemble for every value of `k` from 20 to 40.

	export k
	for k in {20..40}; do
		mkdir k$k
		abyss-pe -C k$k name=ecoli in=../reads.fa
	done
	abyss-fac k*/ecoli-contigs.fa

The default maximum value for `k` is 64. This limit may be changed at
compile time using the `--enable-maxk` option of configure. It may be
decreased to 32 to decrease memory usage or increased to 96.

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

For example, to submit an array of jobs to assemble every odd value of
`k` between 51 and 63 using 64 processes for each job:

	mkdir k{51..63}
	qsub -N ecoli -pe openmpi 64 -t 51-63:2 \
		<<<'abyss-pe -C k$SGE_TASK_ID in=/data/reads.fa'

Assembly Parameters
===================

Parameters of the driver script, `abyss-pe`

 * `a`: maximum number of branches of a bubble [`2`]
 * `b`: maximum length of a bubble (bp) [`10000`]
 * `c`: minimum mean k-mer coverage of a unitig [`sqrt(median)`]
 * `d`: allowable error of a distance estimate (bp) [`6`]
 * `e`: minimum erosion k-mer coverage [`sqrt(median)`]
 * `E`: minimum erosion k-mer coverage per strand [`1`]
 * `j`: number of threads [`2`]
 * `k`: size of k-mer (when `K` is not set) or the span of a k-mer pair (when `K` is set)
 * `K`: the length of a single k-mer in a k-mer pair (bp)
 * `l`: minimum alignment length of a read (bp) [`k`]
 * `m`: minimum overlap of two unitigs (bp) [`30`]
 * `n`: minimum number of pairs required for building contigs [`10`]
 * `N`: minimum number of pairs required for building scaffolds [`n`]
 * `p`: minimum sequence identity of a bubble [`0.9`]
 * `q`: minimum base quality [`3`]
 * `s`: minimum unitig size required for building contigs (bp) [`200`]
 * `S`: minimum contig size required for building scaffolds (bp) [`s`]
 * `t`: minimum tip size (bp) [`2k`]
 * `v`: use `v=-v` for verbose logging, `v=-vv` for extra verbose [`disabled`]

Please see the
[abyss-pe](http://manpages.ubuntu.com/abyss-pe.1.html)
manual page for more information on assembly parameters.

Possibly, `abyss-pe` parameters can have same names as existing environment variables'. The parameters then cannot be used until the environment variables are unset. To detect such occasions, run the command:

	abyss-pe env [options]

Above command will report all `abyss-pe` parameters that are set from various origins. However it will not operate ABySS programs.

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

For a flowchart showing the relationship between these programs,
see doc/flowchart.pdf.

Export to SQLite Database
=========================

ABySS has a built-in support for SQLite database. With this option activated, it exports log values into a SQLite file and/or `.csv` files at runtime.

## Activating the functionality

Download SQLite [here](http://www.sqlite.org/download.html) and install. (See [Quick Start](#quick-start) for details)

To compile ABySS with SQLite, add configure flag `--with-sqlite` to the steps in [Compiling ABySS from GiHub](#compiling-abyss-from-github) / [Compiling ABySS from source](#compiling-abyss-from-source).

	./configure [other options] --with-sqlite=/path/to/sqlite3/
	make
	sudo make install

## Database parameters
Of `abyss-pe`:
 * `db`: path to SQLite repository file [`$(name).sqlite`]
 * `species`: name of species to archive [ ]
 * `strain`: name of strain to archive [ ]
 * `library`: name of library to archive [ ]

For example, to export data of species 'Ecoli', strain 'O121' and library 'pe200' into your SQLite database repository named '/abyss/test.sqlite':

	abyss-pe db=/abyss/test.sqlite species=Ecoli strain=O121 library=pe200 [other options]

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
Steven JM Jones, and İnanç Birol.
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

- **[Shaun Jackman](http://sjackman.ca)**
  — [GitHub/sjackman](https://github.com/sjackman)
  — [@sjackman](https://twitter.com/sjackman)
- **Tony Raymond** — [GitHub/traymond](https://github.com/traymond)
- **Ben Vandervalk** — [GitHub/benvvalk ](https://github.com/benvvalk)
- **Jared Simpson** — [GitHub/jts](https://github.com/jts)

Supervised by [**Dr. İnanç Birol**](http://www.bcgsc.ca/faculty/inanc-birol).

Copyright 2014 Canada's Michael Smith Genome Sciences Centre

[![githalytics.com](https://cruel-carlota.pagodabox.com/af4811df3b40b7d096f6085db2969f0e "githalytics.com")](http://githalytics.com/sjackman/abyss)
