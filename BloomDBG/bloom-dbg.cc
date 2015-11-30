#include "config.h"

#include "BloomDBG/bloom-dbg.h"
#include "BloomDBG/HashAgnosticCascadingBloom.h"
#include "Common/Kmer.h"
#include "Common/StringUtil.h"
#include "Common/Options.h"
#include "lib/bloomfilter-9061f087d8714109b865415f2ac05e4796d0cd74/BloomFilter.hpp"

#include <getopt.h>
#include <iostream>
#include <sstream>
#include <cstdlib>

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
"Usage: " PROGRAM " -b <bloom_size> -H <bloom_hashes> -k <kmer_size> [options] \\\n"
"    <FASTQ> [FASTQ]... > assembly.fasta\n"
"\n"
"Perform a de Bruijn graph assembly of the given FASTQ files.\n"
"\n"
"Options:\n"
"\n"
"  -b  --bloom-size=N         Bloom filter memory size with suffix\n"
"                             of 'k', 'M', or 'G' [required]\n"
"  -H  --num-hashes=N         number of Bloom filter hash functions\n"
"                             [required]\n"
"      --help                 display this help and exit\n"
"  -j, --threads=N            use N parallel threads [1]\n"
"  -k, --kmer=N               the size of a k-mer [required]\n"
"  -v, --verbose              display verbose output\n"
"      --version              output version information and exit\n"
"\n"
"Example:\n"
"\n"
"  Assemble a 100 Mbp genome with 50X coverage and a k-mer size of 50bp:\n"
"\n"
"  $ " PROGRAM " -C50 -G100M -k50 reads1.fq.gz reads2.fq.gz > assembly.fa\n"
"\n"
"Report bugs to <" PACKAGE_BUGREPORT ">.\n";

namespace opt {

	/** Bloom filter size (in bits) */
	static size_t bloomSize = 0;

	/** number of cascading Bloom filter levels */
	static unsigned bloomLevels = 2;

	/** num Bloom filter hash functions */
	static unsigned numHashes = 0;

	/** The number of parallel threads. */
	static unsigned threads = 1;

	/** The size of a k-mer. */
	static unsigned k;

}

static const char shortopts[] = "b:H:j:k:v";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
	{ "help",             no_argument, NULL, OPT_HELP },
	{ "threads",          required_argument, NULL, 'j' },
	{ "kmer",             required_argument, NULL, 'k' },
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
			opt::bloomSize = SIToBytes(arg); break;
		  case 'H':
			arg >> opt::numHashes; break;
		  case 'j':
			arg >> opt::threads; break;
		  case 'k':
			arg >> opt::k; break;
		  case 'v':
			++opt::verbose; break;
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

	if (opt::bloomSize == 0) {
		cerr << PROGRAM ": missing mandatory option `-b'\n";
		die = true;
	}

	if (opt::numHashes == 0) {
		cerr << PROGRAM ": missing mandatory option `-H'\n";
		die = true;
	}

	if (opt::k == 0) {
		cerr << PROGRAM ": missing mandatory option `-k'\n";
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
	if (opt::threads > 0)
		omp_set_num_threads(opt::threads);
#endif

	/* set global variable for k-mer length */
	Kmer::setLength(opt::k);

	/* use cascading Bloom filter to remove error k-mers */
	const size_t bitsPerByte = 8;
	HashAgnosticCascadingBloom cascadingBloom(
		opt::bloomSize * bitsPerByte / opt::bloomLevels,
		opt::numHashes, opt::bloomLevels, opt::k);

	/* load reads into Bloom filter */
	for (; optind < argc; ++optind) {
		BloomDBG::loadFile(cascadingBloom, opt::numHashes,
			opt::k, argv[optind], opt::verbose);
	}

	return EXIT_SUCCESS;
}
