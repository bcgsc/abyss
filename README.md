[![Release](https://img.shields.io/github/release/bcgsc/abyss.svg)](https://github.com/bcgsc/abyss/releases)
[![Downloads](https://img.shields.io/github/downloads/bcgsc/abyss/total?logo=github)](https://github.com/bcgsc/abyss/releases/download/2.2.3/abyss-2.2.3.tar.gz)
[![Conda](https://img.shields.io/conda/dn/bioconda/abyss?label=Conda)](https://anaconda.org/bioconda/abyss)
[![Issues](https://img.shields.io/github/issues/bcgsc/abyss.svg)](https://github.com/bcgsc/abyss/issues)

ABySS
================================================================================

ABySS is a *de novo* sequence assembler intended for short paired-end reads and genomes of all sizes.

Please [cite our papers](#citation).

News
================================================================================
3 May 2019

Looking for a fun & worthy challenge? Think you can contribute code to this project? 
Join our team of developers! We are currently looking for C++ bioinformatics programmers.
Inquire for staff, graduate student, and postdoctoral positions. 
[Contact the project lead (Inanc Birol)](mailto:ibirol@bcgsc.ca?Subject=ABySS%20developer%20position) 

Contents
================================================================================

* [Installation](#installation)
	* [Install ABySS using Conda](#install-abyss-using-conda)
	* [Install ABySS using Homebrew](#install-abyss-using-homebrew)
	* [Install ABySS using apt-get](#install-abyss-using-apt-get)
	* [Install ABySS on Windows](#install-abyss-on-windows)
* [Dependencies](#dependencies)
	* [Dependencies for linked reads](#dependencies-for-linked-reads)
	* [Optional dependencies](#optional-dependencies)
* [Compiling ABySS from source](#compiling-abyss-from-source)
* [Before starting an assembly](#before-starting-an-assembly)
* [Modes](#modes)
	* [Bloom filter mode](#bloom-filter-mode)
	* [MPI mode (legacy)](#mpi-mode-legacy)
* [Examples](#examples)
	* [Assemble a small synthetic data set](#assemble-a-small-synthetic-data-set)
	* [Assembling a paired-end library](#assembling-a-paired-end-library)
	* [Assembling multiple libraries](#assembling-multiple-libraries)
	* [Scaffolding](#scaffolding)
	* [Scaffolding with linked reads](#scaffolding-with-linked-reads)
	* [Rescaffolding with long sequences](#rescaffolding-with-long-sequences)
	* [Assembling using a paired de Bruijn graph](#assembling-using-a-paired-de-bruijn-graph)
	* [Assembling a strand-specific RNA-Seq library](#assembling-a-strand-specific-rna-seq-library)
* [Optimizing the parameters k and kc](#optimizing-the-parameters-k-and-kc)
* [Running ABySS on a cluster](#running-abyss-on-a-cluster)
* [Using the DIDA alignment framework](#using-the-dida-alignment-framework)
* [Assembly Parameters](#assembly-parameters)
* [ABySS programs](#abyss-programs)
* [Export to SQLite Database](#export-to-sqlite-database)
	* [Database parameters](#database-parameters)
	* [Helper programs](#helper-programs)
* [Citation](#citation)
* [Related Publications](#related-publications)
* [Support](#support)
* [Authors](#authors)

Installation
================================================================================

## Install ABySS using Conda (recommended)

If you have the [Conda](https://docs.conda.io/en/latest/) package manager (Linux, MacOS) installed, run:

	conda install -c bioconda abyss

Or you can install ABySS in a dedicated environment:

    conda create -n abyss-env
    conda activate abyss-env
    conda install -c bioconda abyss

## Install ABySS using Homebrew

If you have the [Homebrew](https://brew.sh) package manager (Linux, MacOS) installed, run:

	brew install abyss

## Install ABySS on Windows

Install [Windows Subsystem for Linux](https://docs.microsoft.com/en-us/windows/wsl/) from which you can run Conda or Homebrew installation.

Dependencies
============

## Dependencies for linked reads

- [ARCS](https://github.com/bcgsc/arcs) for scaffolding.
- [Tigmint](https://github.com/bcgsc/tigmint) for correcting assembly errors.

These can be installed through Conda:

	conda install -c bioconda arcs tigmint

Or Homebrew:

	brew install brewsci/bio/arcs brewsci/bio/links-scaffolder

## Optional dependencies

- [pigz](https://zlib.net/pigz/) for parallel gzip.
- [samtools](https://samtools.github.io) for reading BAM files.
- [zsh](https://sourceforge.net/projects/zsh/) for reporting time and memory usage.

Conda:

	conda install -c bioconda samtools
	conda install -c conda-forge pigz zsh

Homebrew:

	brew install pigz samtools zsh

Compiling ABySS from source
================================================================================

When compiling ABySS from source the following tools are
required:

* [Autoconf](http://www.gnu.org/software/autoconf)
* [Automake](http://www.gnu.org/software/automake)

ABySS requires a C++ compiler that supports
[OpenMP](http://www.openmp.org) such as [GCC](http://gcc.gnu.org).

The following libraries are required:

* [Boost](http://www.boost.org/)
* [Open MPI](http://www.open-mpi.org)
* [sparsehash](https://code.google.com/p/sparsehash/)

Conda:

	conda install -c conda-forge boost openmpi
	conda install -c bioconda google-sparsehash

It is also helpful to install the compilers Conda package that automatically passes the correct compiler flags to use the available Conda packages:

	conda install -c conda-forge compilers

Homebrew:

	brew install boost open-mpi google-sparsehash

ABySS will receive an error when compiling with Boost 1.51.0 or 1.52.0
since they contain a bug. Later versions of Boost compile without error.

To compile, run the following:

	./autogen.sh
	mkdir build
	cd build
	../configure --prefix=/path/to/abyss
	make
	make install

You may also pass the following flags to `configure` script:

	--with-boost=PATH
	--with-mpi=PATH
	--with-sqlite=PATH
	--with-sparsehash=PATH

Where PATH is the path to the directory containing the corresponding dependencies. This should only be necessary if `configure` doesn't find the dependencies by default. If you are using Conda, PATH would be the path to the Conda installation. SQLite and MPI are optional dependencies.

The above steps install ABySS at the provided path, in this case `/path/to/abyss`.
Not specifying `--prefix` would install in `/usr/local`, which requires
sudo privileges when running `make install`.

ABySS requires a modern compiler such as GCC 6 or greater. If you have multiple
versions of GCC installed, you can specify a different compiler:

	../configure CC=gcc-10 CXX=g++-10

While OpenMPI is assumed by default you can switch to LAM/MPI or MPICH
using:

        ../configure --enable-lammpi
        ../configure --enable-mpich

The default maximum k-mer size is 192 and may be decreased to reduce
memory usage or increased at compile time. This value must be a
multiple of 32 (i.e. 32, 64, 96, 128, etc):

	../configure --enable-maxk=160

If you encounter compiler warnings that are not critical, you can allow the compilation to continue:

	../configure --disable-werror

To run ABySS, its executables should be found in your `PATH` environment variable. If you
installed ABySS in `/opt/abyss`, add `/opt/abyss/bin` to your `PATH`:

	PATH=/opt/abyss/bin:$PATH

Before starting an assembly
================================================================================

ABySS stores temporary files in `TMPDIR`, which is `/tmp` by default on most systems. If your default temporary disk volume is too small, set `TMPDIR` to a larger volume, such as `/var/tmp` or your home directory.

	export TMPDIR=/var/tmp

Modes
================================================================================

## Bloom filter mode

The recommended mode of running ABySS is the Bloom filter mode. Specifying
the Bloom filter memory budget with the `B` parameter enables this mode, which can
reduce memory consumption by ten-fold compared to the MPI mode. `B` may be specified
with unit suffixes 'k' (kilobytes), 'M' (megabytes), 'G' (gigabytes). If no units
are specified bytes are assumed. Internally, the Bloom filter assembler allocates
the entire memory budget (`B * 8/9`) to a Counting Bloom filter, and an additional
(`B/9`) memory to another Bloom filter that is used to track k-mers that have previously
been included in contigs.

A good value for `B` depends on a number of factors, but primarily on the
genome being assembled. A general guideline is:

P. glauca (~20Gbp): `B=500G`
H. sapiens (~3.1Gbp): `B=50G`
C. elegans (~101Mbp): `B=2G`

For other genome sizes, the value for `B` can be interpolated. Note that
there is no downside to using larger than necessary `B` value, except for
the memory required. To make sure you have selected a correct `B` value,
inspect the standard error log of the assembly process and ensure that the
reported FPR value under `Counting Bloom filter stats` is 5% or less. This
requires using verbosity level 1 with `v=-v` option.

## MPI mode (legacy)

This mode is legacy and we do not recommend running ABySS with it.
To run ABySS in the MPI mode, you need to specify the `np` parameter,
which specifies the number of processes to use for the parallel MPI job.
Without any MPI configuration, this will allow you to use multiple cores
on a single machine. To use multiple machines for assembly, you must create
a `hostfile` for `mpirun`, which is described in the `mpirun` man page.

*Do not* run `mpirun -np 8 abyss-pe`. To run ABySS with 8 threads, use
`abyss-pe np=8`. The `abyss-pe` driver script will start the MPI
process, like so: `mpirun -np 8 ABYSS-P`.

The paired-end assembly stage is multithreaded, but must run on a
single machine. The number of threads to use may be specified with the
parameter `j`. The default value for `j` is the value of `np`.

Examples
================================================================================

## Assemble a small synthetic data set

	wget http://www.bcgsc.ca/platform/bioinfo/software/abyss/releases/1.3.4/test-data.tar.gz
	tar xzvf test-data.tar.gz
	abyss-pe k=25 name=test B=1G \
		in='test-data/reads1.fastq test-data/reads2.fastq'

Calculate assembly contiguity statistics:

	abyss-fac test-unitigs.fa

## Assembling a paired-end library

To assemble paired reads in two files named `reads1.fa` and
`reads2.fa` into contigs in a file named `ecoli-contigs.fa`, run the
command:

	abyss-pe name=ecoli k=96 B=2G in='reads1.fa reads2.fa'

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

## Assembling multiple libraries

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

	abyss-pe k=96 B=2G name=ecoli lib='pea peb' \
		pea='pea_1.fa pea_2.fa' peb='peb_1.fa peb_2.fa' \
		se='se1.fa se2.fa'

The empirical distribution of fragment sizes will be stored in two
files named `pea-3.hist` and `peb-3.hist`. These files may be
plotted to check that the empirical distribution agrees with the
expected distribution. The assembled contigs will be stored in
`${name}-contigs.fa`.

## Scaffolding

Long-distance mate-pair libraries may be used to scaffold an assembly.
Specify the names of the mate-pair libraries using the parameter `mp`.
The scaffolds will be stored in the file `${name}-scaffolds.fa`.
Here's an example of assembling a data set with two paired-end
libraries and two mate-pair libraries. Note that the names of the libraries
(`pea`, `peb`, `mpa`, `mpb`) are arbitrary.

	abyss-pe k=96 B=2G name=ecoli lib='pea peb' mp='mpc mpd' \
		pea='pea_1.fa pea_2.fa' peb='peb_1.fa peb_2.fa' \
		mpc='mpc_1.fa mpc_2.fa' mpd='mpd_1.fa mpd_2.fa'

The mate-pair libraries are used only for scaffolding and do not
contribute towards the consensus sequence.

## Scaffolding with linked reads

ABySS can scaffold using linked reads from 10x Genomics Chromium. The barcodes must first be extracted from the read sequences and added to the `BX:Z` tag of the FASTQ header, typically using the `longranger basic` command of [Long Ranger](https://support.10xgenomics.com/genome-exome/software/overview/welcome) or [EMA preproc](https://github.com/arshajii/ema#readme). The linked reads are used to correct assembly errors, which requires that [Tigmint](https://github.com/bcgsc/tigmint). The linked reads are also used for scaffolding, which requires [ARCS](https://github.com/bcgsc/arcs). See [Dependencies](#dependencies) for installation instructions.

ABySS can combine paired-end, mate-pair, and linked-read libraries. The `pe` and `lr` libraries will be used to build the de Bruijn graph. The `mp` libraries will be used for paired-end/mate-pair scaffolding. The `lr` libraries will be used for misassembly correction using Tigmint and scaffolding using ARCS.

	abyss-pe k=96 B=2G name=hsapiens \
		pe='pea' pea='lra.fastq.gz' \
		mp='mpa' mpa='lra.fastq.gz' \
		lr='lra' lra='lra.fastq.gz'

ABySS performs better with a mixture of paired-end, mate-pair, and linked reads, but it is possible to assemble only linked reads using ABySS, though this mode of operation is experimental.

	abyss-pe k=96 name=hsapiens lr='lra' lra='lra.fastq.gz'

## Rescaffolding with long sequences

Long sequences such as RNA-Seq contigs can be used to rescaffold an
assembly. Sequences are aligned using BWA-MEM to the assembled
scaffolds. Additional scaffolds are then formed between scaffolds that
can be linked unambiguously when considering all BWA-MEM alignments.

Similar to scaffolding, the names of the datasets can be specified with
the `long` parameter. These scaffolds will be stored in the file
`${name}-long-scaffs.fa`. The following is an example of an assembly with PET, MPET and an RNA-Seq assembly. Note that the names of the libraries are arbitrary.

	abyss-pe k=96 B=2G name=ecoli lib='pe1 pe2' mp='mp1 mp2' long='longa' \
		pe1='pe1_1.fa pe1_2.fa' pe2='pe2_1.fa pe2_2.fa' \
		mp1='mp1_1.fa mp1_2.fa' mp2='mp2_1.fa mp2_2.fa' \
		longa='longa.fa'

## Assembling using a paired de Bruijn graph

Assemblies may be performed using a _paired de Bruijn graph_ instead
of a standard de Bruijn graph. In paired de Bruijn graph mode, ABySS
uses _k-mer pairs_ in place of k-mers, where each k-mer pair consists of
two equal-size k-mers separated by a fixed distance. A k-mer pair
is functionally similar to a large k-mer spanning the breadth of the k-mer
pair, but uses less memory because the sequence in the gap is not stored.
To assemble using paired de Bruijn graph mode, specify both individual
k-mer size (`K`) and k-mer pair span (`k`). For example, to assemble E.
coli with a individual k-mer size of 16 and a k-mer pair span of 96:

	abyss-pe name=ecoli K=16 k=96 in='reads1.fa reads2.fa'

In this example, the size of the intervening gap between k-mer pairs is
64 bp (96 - 2\*16). Note that the `k` parameter takes on a new meaning
in paired de Bruijn graph mode. `k` indicates kmer pair span in
paired de Bruijn graph mode (when `K` is set), whereas `k` indicates
k-mer size in standard de Bruijn graph mode (when `K` is not set).

## Assembling a strand-specific RNA-Seq library

Strand-specific RNA-Seq libraries can be assembled such that the
resulting unitigs, contigs and scaffolds are oriented correctly with
respect to the original transcripts that were sequenced. In order to
run ABySS in strand-specific mode, the `SS` parameter must be used as
in the following example:

	abyss-pe name=SS-RNA B=2G k=96 in='reads1.fa reads2.fa' SS=--SS

The expected orientation for the read sequences with respect to the
original RNA is RF. i.e. the first read in a read pair is always in
reverse orientation.

Optimizing the parameters k and kc
================================================================================

It is standard practice when running ABySS to run multiple assemblies
to find the optimal values for the `k` and `kc` parameters. `k` determines
the k-mer size in the de Bruijn Graph, and `kc` is the k-mer minimum coverage
multiplicity cutoff, which filters out erroneous k-mers. The range in which `k`
should be tested depends on the read size and read coverage.

A rough indicator is, for 2x150bp reads and 40x coverage, the right `k` value is often around 70 to 90. For 2x250bp reads and 40x coverage, the right value might be around 110 to 140.

For `kc`, 2 is most often a good value, but can go as high as 4.

The following shell snippet will assemble for `k` values 2 and 3, and every eighth value of `k` from 50 to 90. In the end, we calculate the contiguity statistics, as a proxy for identifying the optimal assembly. Other metrics can be used, as needed.

	for kc in 2 3; do
		for k in `seq 50 8 90`; do
			mkdir k${k}-kc${kc}
			abyss-pe -C k${k}-kc${kc} name=ecoli B=2G k=$k kc=$kc in=../reads.fa
		done
	done
	abyss-fac k*/ecoli-scaffolds.fa

The default maximum value for `k` is 192. This limit may be changed at
compile time using the `--enable-maxk` option of configure. It may be
decreased to 32 to decrease memory usage or increased to larger values.

Running ABySS on a cluster
================================================================================

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
================================================================================

ABySS supports the use of DIDA (Distributed Indexing Dispatched Alignment),
an MPI-based framework for computing sequence alignments in parallel across
multiple machines. The DIDA software must be separately downloaded and
installed from http://www.bcgsc.ca/platform/bioinfo/software/dida. In
comparison to the standard ABySS alignment stages which are constrained
to a single machine, DIDA offers improved performance and the ability to
scale to larger targets. Please see the DIDA section of the abyss-pe man
page (in the `doc` subdirectory) for details on usage.

Assembly Parameters
================================================================================

Parameters of the driver script, `abyss-pe`

 * `a`: maximum number of branches of a bubble [`2`]
 * `b`: maximum length of a bubble (bp) [`""`]
 * `B`: Bloom filter size (e.g. "100M")
 * `c`: minimum mean k-mer coverage of a unitig [`sqrt(median)`]
 * `d`: allowable error of a distance estimate (bp) [`6`]
 * `e`: minimum erosion k-mer coverage [`round(sqrt(median))`]
 * `E`: minimum erosion k-mer coverage per strand [1 if `sqrt(median) > 2` else 0]
 * `G`: genome size, used to calculate NG50
 * `H`: number of Bloom filter hash functions [`4`]
 * `j`: number of threads [`2`]
 * `k`: size of k-mer (when `K` is not set) or the span of a k-mer pair (when `K` is set)
 * `kc`: minimum k-mer count threshold for Bloom filter assembly [`2`]
 * `K`: the length of a single k-mer in a k-mer pair (bp)
 * `l`: minimum alignment length of a read (bp) [`40`]
 * `m`: minimum overlap of two unitigs (bp) [`0` (interpreted as `k - 1`) if `mp` is provided or if `k<=50`, otherwise `50`]
 * `n`: minimum number of pairs required for building contigs [`10`]
 * `N`: minimum number of pairs required for building scaffolds [`15-20`]
 * `np`: number of MPI processes [`1`]
 * `p`: minimum sequence identity of a bubble [`0.9`]
 * `q`: minimum base quality [`3`]
 * `s`: minimum unitig size required for building contigs (bp) [`1000`]
 * `S`: minimum contig size required for building scaffolds (bp) [`100-5000`]
 * `t`: maximum length of blunt contigs to trim [`k`]
 * `v`: use `v=-v` for verbose logging, `v=-vv` for extra verbose
 * `x`: spaced seed (Bloom filter assembly only)
 * `lr_s`: minimum contig size required for building scaffolds with linked reads (bp) [`S`]
 * `lr_n`: minimum number of barcodes required for building scaffolds with linked reads [`10`]

Environment variables
================================================================================

`abyss-pe` configuration variables may be set on the command line or from the environment, for example with `export k=96`. It can happen that `abyss-pe` picks up such variables from your environment that you had not intended, and that can cause trouble. To troubleshoot that situation, use the `abyss-pe env` command to print the values of all the `abyss-pe` configuration variables:

	abyss-pe env [options]

ABySS programs
================================================================================

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
 * `abyss-rresolver`: resolve repeats using short reads

This [flowchart](https://github.com/bcgsc/abyss/blob/master/doc/flowchart.pdf) shows the ABySS assembly pipeline and its intermediate files.

Export to SQLite Database
================================================================================

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

Citation
================================================================================

## [ABySS 2.0](http://doi.org/10.1101/gr.214346.116)

Shaun D Jackman, Benjamin P Vandervalk, Hamid Mohamadi, Justin Chu, Sarah Yeo, S Austin Hammond, Golnaz Jahesh, Hamza Khan, Lauren Coombe, René L Warren, and Inanc Birol (2017).
**ABySS 2.0: Resource-efficient assembly of large genomes using a Bloom filter**.
*Genome research*, 27(5), 768-777.
[doi:10.1101/gr.214346.116](http://doi.org/10.1101/gr.214346.116)

## [ABySS](http://genome.cshlp.org/content/19/6/1117)

Simpson, Jared T., Kim Wong, Shaun D. Jackman, Jacqueline E. Schein,
Steven JM Jones, and Inanc Birol (2009).
**ABySS: a parallel assembler for short read sequence data**.
*Genome research*, 19(6), 1117-1123.
[doi:10.1101/gr.089532.108](http://dx.doi.org/10.1101/gr.089532.108)

Related Publications
================================================================================

## [RResolver](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-022-04790-z)

Vladimir Nikolić, Amirhossein Afshinfard, Justin Chu, Johnathan Wong, Lauren Coombe, Ka Ming Nip, René L. Warren & Inanç Birol (2022).
**RResolver: efficient short-read repeat resolution within ABySS**.
*BMC Bioinformatics* 23, Article number: 246 (2022).
[doi:10.1186/s12859-022-04790-z](https://doi.org/10.1186/s12859-022-04790-z)

## [Trans-ABySS](http://www.nature.com/nmeth/journal/v7/n11/abs/nmeth.1517.html)

Robertson, Gordon, Jacqueline Schein, Readman Chiu, Richard Corbett,
Matthew Field, Shaun D. Jackman, Karen Mungall, et al (2010).
**De novo assembly and analysis of RNA-seq data**.
*Nature methods*, 7(11), 909-912.
[doi:10.1038/10.1038/nmeth.1517](http://dx.doi.org/10.1038/nmeth.1517)

## [ABySS-Explorer](http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=5290690)

Nielsen, Cydney B., Shaun D. Jackman, Inanc Birol, and Steven JM Jones (2009).
**ABySS-Explorer: visualizing genome sequence assemblies**.
*IEEE Transactions on Visualization and Computer Graphics*, 15(6), 881-888.
[doi:10.1109/TVCG.2009.116](http://dx.doi.org/10.1109/TVCG.2009.116)

Support
================================================================================

[Create a new issue on GitHub.](https://github.com/bcgsc/abyss/issues)

[Ask a question on Biostars.](https://www.biostars.org/tag/abyss/)

Subscribe to the [ABySS mailing list](http://groups.google.com/group/abyss-users), <abyss-users@googlegroups.com>.

For questions related to transcriptome assembly, contact the [Trans-ABySS mailing list](http://groups.google.com/group/trans-abyss), <trans-abyss@googlegroups.com>.

Authors
================================================================================

+ **[Shaun Jackman](http://sjackman.ca)** - [GitHub/sjackman](https://github.com/sjackman) - [@sjackman](https://twitter.com/sjackman)
+ **Tony Raymond** - [GitHub/traymond](https://github.com/traymond)
+ **Ben Vandervalk** - [GitHub/benvvalk ](https://github.com/benvvalk)
+ **Jared Simpson** - [GitHub/jts](https://github.com/jts)
+ **Johnathan Wong** - [GitHub/jowong4](https://github.com/jowong4)
+ **Vladimir Nikolić** - [GitHub/vlad0x00](https://github.com/vlad0x00)

Supervised by [**Dr. Inanc Birol**](http://www.bcgsc.ca/faculty/inanc-birol).

Copyright 2016 Canada's Michael Smith Genome Sciences Centre
