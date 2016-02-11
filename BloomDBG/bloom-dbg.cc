#include "config.h"

#include "BloomDBG/bloom-dbg.h"
#include "BloomDBG/HashAgnosticCascadingBloom.h"
#include "Common/Kmer.h"
#include "Common/StringUtil.h"
#include "Common/Options.h"
#include "DataLayer/Options.h"
#include "lib/bloomfilter-2dfba08d120d7659e8c75cf5c501b3b9040e98cb/BloomFilter.hpp"

#include <getopt.h>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <iomanip>
#include <cstring>

#if _OPENMP
# include <omp.h>
#endif

using namespace std;

#define PROGRAM "abyss-bloom-dbg"

static const char VERSION_MESSAGE[] =
	PROGRAM " (" PACKAGE_NAME ") " VERSION "\n"
	"Written by Ben Vandervalk, Shaun Jackman, Hamid Mohamadi,\n"
	"Justin Chu, and Anthony Raymond.\n"
	"\n"
	"Copyright 2015 Canada's Michael Smith Genome Science Centre\n";

static const char USAGE_MESSAGE[] =
"Usage: " PROGRAM " -b <bloom_size> -H <bloom_hashes> -k <kmer_size> \\\n"
"    -G <genome_size> [options] <FASTQ> [FASTQ]... > assembly.fasta\n"
"\n"
"Perform a de Bruijn graph assembly of the given FASTQ files.\n"
"\n"
"Options:\n"
"\n"
"  -b  --bloom-size=N         Bloom filter memory size with unit suffix\n"
"                             'k', 'M', or 'G' [required]\n"
"  -c, --min-coverage=N       kmer coverage threshold for error correction [2]\n"
"      --chastity             discard unchaste reads [default]\n"
"      --no-chastity          do not discard unchaste reads\n"
"  -g  --graph=FILE           write de Bruijn graph to FILE (GraphViz)\n"
"  -G  --genome-size=N        approx genome size with unit suffix\n"
"                             'k', 'M', or 'G' [required]\n"
"  -H  --num-hashes=N         number of Bloom filter hash functions\n"
"                             [required]\n"
"      --help                 display this help and exit\n"
"  -j, --threads=N            use N parallel threads [1]\n"
"      --trim-masked          trim masked bases from the ends of reads\n"
"      --no-trim-masked       do not trim masked bases from the ends\n"
"                             of reads [default]\n"
"  -k, --kmer=N               the size of a k-mer [required]\n"
"  -q, --trim-quality=N       trim bases from the ends of reads whose\n"
"                             quality is less than the threshold\n"
"  -Q, --mask-quality=N       mask all low quality bases as `N'\n"
"      --standard-quality     zero quality is `!' (33), typically\n"
"                             for FASTQ and SAM files [default]\n"
"      --illumina-quality     zero quality is `@' (64), typically\n"
"                             for qseq and export files\n"
"  -s, --spaced-seed=STR      bitmask indicating k-mer positions to be\n"
"                             ignored during hashing [default is string\n"
"                             of '1's]\n"
"  -v, --verbose              display verbose output\n"
"      --version              output version information and exit\n"
"\n"
"Example:\n"
"\n"
"  Assemble a 100 Mbp genome using a k-mer size of 50bp. Allocate a 1GB\n"
"  Bloom filter with 2 hash functions and require that a k-mer\n"
"  occurs 3 times or more to be included in the assembly. (The k-mer\n"
"  count threshold filters out k-mers due to sequencing errors.)\n"
"\n"
"  $ " PROGRAM " -G100M -k50 -b1G -H2 -c3 reads1.fq.gz reads2.fq.gz > assembly.fa\n"
"\n"
"Report bugs to <" PACKAGE_BUGREPORT ">.\n";

/** Assembly params (stores command-line options) */
BloomDBG::AssemblyParams params;

static const char shortopts[] = "b:c:g:G:H:j:k:q:Q:s:v";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
	{ "bloom-size",       required_argument, NULL, 'b' },
	{ "min-coverage",     required_argument, NULL, 'c' },
	{ "chastity",         no_argument, &opt::chastityFilter, 1 },
	{ "no-chastity",      no_argument, &opt::chastityFilter, 0 },
	{ "graph",            required_argument, NULL, 'g' },
	{ "genome-size",      required_argument, NULL, 'G' },
	{ "num-hashes",       required_argument, NULL, 'H' },
	{ "help",             no_argument, NULL, OPT_HELP },
	{ "threads",          required_argument, NULL, 'j' },
	{ "trim-masked",      no_argument, &opt::trimMasked, 1 },
	{ "no-trim-masked",   no_argument, &opt::trimMasked, 0 },
	{ "kmer",             required_argument, NULL, 'k' },
	{ "trim-quality",     required_argument, NULL, 'q' },
	{ "mask-quality",     required_argument, NULL, 'Q' },
	{ "standard-quality", no_argument, &opt::qualityOffset, 33 },
	{ "illumina-quality", no_argument, &opt::qualityOffset, 64 },
	{ "spaced-seed",      no_argument, NULL, 's' },
	{ "verbose",          no_argument, NULL, 'v' },
	{ "version",          no_argument, NULL, OPT_VERSION },
	{ NULL, 0, NULL, 0 }
};

/**
 * Create a de novo genome assembly using a Bloom filter de
 * Bruijn graph.
 */
int main(int argc, char** argv)
{
	bool die = false;

	for (int c; (c = getopt_long(argc, argv,
					shortopts, longopts, NULL)) != -1;) {
		istringstream arg(optarg != NULL ? optarg : "");
		switch (c) {
		  case '?':
			die = true; break;
		  case 'b':
			params.bloomSize = SIToBytes(arg); break;
		  case 'c':
			arg >> params.minCov; break;
		  case 'g':
			arg >> params.graphPath; break;
		  case 'G':
			  params.genomeSize = fromSI(arg); break;
		  case 'H':
			arg >> params.numHashes; break;
		  case 'j':
			arg >> params.threads; break;
		  case 'k':
			arg >> params.k; break;
		  case 'q':
			arg >> opt::qualityThreshold; break;
		  case 's':
			arg >> params.spacedSeed; break;
		  case 'Q':
			arg >> opt::internalQThreshold; break;
		  case 'v':
			++params.verbose; break;
		  case OPT_HELP:
			cout << USAGE_MESSAGE;
			exit(EXIT_SUCCESS);
		  case OPT_VERSION:
			cout << VERSION_MESSAGE;
			exit(EXIT_SUCCESS);
		}
		if (optarg != NULL && (!arg.eof() || arg.fail())) {
			cerr << PROGRAM ": invalid option: `-"
				<< (char)c << optarg << "'\n";
			exit(EXIT_FAILURE);
		}
	}

	if (params.bloomSize == 0) {
		cerr << PROGRAM ": missing mandatory option `-b'\n";
		die = true;
	}

	if (params.genomeSize == 0) {
		cerr << PROGRAM ": missing mandatory option `-G'\n";
		die = true;
	}

	if (params.numHashes == 0) {
		cerr << PROGRAM ": missing mandatory option `-H'\n";
		die = true;
	}

	if (params.k == 0) {
		cerr << PROGRAM ": missing mandatory option `-k'\n";
		die = true;
	}

	if (params.spacedSeed.empty()) {
		/* spaced seed defaults to all '1's */
		params.spacedSeed.reserve(params.k);
		for (size_t i = 0; i < params.k; ++i)
			params.spacedSeed.push_back('1');
	} else if (params.spacedSeed.length() != params.k) {
		cerr << PROGRAM ": spaced seed must be exactly k bits long\n";
		die = true;
	} else if (params.spacedSeed.find_first_not_of("01") != string::npos) {
		cerr << PROGRAM ": spaced seed must contain only '0's or '1's\n";
		die = true;
	}

	if (argc - optind < 1) {
		cerr << PROGRAM ": missing input file arguments\n";
		die = true;
	}

	if (die) {
		cerr << "Try `" << PROGRAM
			<< " --help' for more information.\n";
		exit(EXIT_FAILURE);
	}

#if _OPENMP
	if (params.threads > 0)
		omp_set_num_threads(params.threads);
#endif

	/* set global variable for k-mer length */
	Kmer::setLength(params.k);

	/* BloomFilter class requires size to be a multiple of 64 */
	const size_t bitsPerByte = 8;
	size_t bloomLevelSize = BloomDBG::roundUpToMultiple(
		params.bloomSize * bitsPerByte / params.minCov, (size_t)64);

	/* use cascading Bloom filter to remove error k-mers */
	HashAgnosticCascadingBloom cascadingBloom(
		bloomLevelSize, params.numHashes, params.minCov, params.k);

	/* load reads into Bloom filter */
	for (int i = optind; i < argc; ++i) {
		/*
		 * Debugging feature: If there is a ':'
		 * separating the list of input read files into
		 * two parts, use the first set of files
		 * to load the Bloom filter and the second
		 * set of files for the assembly (read extension).
		 */
		if (strcmp(argv[i],":") == 0) {
			optind = i + 1;
			break;
		}
		BloomDBG::loadFile(cascadingBloom, argv[i], params.spacedSeed,
			params.verbose);
	}
	if (params.verbose)
		cerr << "Bloom filter FPR: " << setprecision(3)
			<< cascadingBloom.FPR() * 100 << "%" << endl;

	/* second pass through FASTA files for assembling */
	BloomDBG::assemble(argc - optind, argv + optind,
		cascadingBloom, params, cout);

	/* generate de Bruijn graph in GraphViz format (optional) */
	if (!params.graphPath.empty()) {
		ofstream graphOut(params.graphPath.c_str());
		assert_good(graphOut, params.graphPath);
		BloomDBG::outputGraph(argc - optind, argv + optind,
			cascadingBloom, params, graphOut);
		assert_good(graphOut, params.graphPath);
		graphOut.close();
		assert_good(graphOut, params.graphPath);
	}

	return EXIT_SUCCESS;
}
