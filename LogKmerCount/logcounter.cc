/**
 * logarithmic k-mer count based on a Bloom filter
 * Copyright 2014 Justin Chu
 */

#include "config.h"

#include <cassert>
#include <getopt.h>
#include <iostream>
#if _OPENMP
# include <omp.h>
#endif

using namespace std;

#define PROGRAM "logcounter"

static const char VERSION_MESSAGE[] =
PROGRAM " (" PACKAGE_NAME ") " VERSION "\n"
"Written by ?.\n"
"\n"
"Copyright 2014 Canada's Michael Smith Genome Science Centre\n";

static const char USAGE_MESSAGE[] =
"Usage: " PROGRAM " [OPTION]... [READS]...\n"
"  -j, --threads=N            use N parallel threads [1]\n"
"  -k, --kmer=N               the size of a k-mer\n"
//"  -b, --bloom-size=N         size of bloom filter [500M]\n"
"      --chastity             discard unchaste reads [default]\n"
"      --no-chastity          do not discard unchaste reads\n"
"      --trim-masked          trim masked bases from the ends of reads\n"
"      --no-trim-masked       do not trim masked bases from the ends\n"
"                             of reads [default]\n"
"  -q, --trim-quality=N       trim bases from the ends of reads whose\n"
"                             quality is less than the threshold\n"
"      --standard-quality     zero quality is `!' (33)\n"
"                             default for FASTQ and SAM files\n"
"      --illumina-quality     zero quality is `@' (64)\n"
"                             default for qseq and export files\n"
"  -v, --verbose              display verbose output\n"
"      --help                 display this help and exit\n"
"      --version              output version information and exit\n"
"\n"
"Report bugs to <" PACKAGE_BUGREPORT ">.\n";

const unsigned g_progressStep = 1000;

namespace opt {
	/** The number of parallel threads. */
	static unsigned threads = 1;

//	/** The size of the bloom filter in bytes. */
//	size_t bloomSize = 500 * 1024 * 1024;

	/** The size of a k-mer. */
	unsigned k;

//	/** Prefix for output files */
//	static string outputPrefix;
}

/** Counters */
//static struct {
//
//} g_count;

static const char shortopts[] = "j:k:q:v";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
//	{ "bloom-size",       required_argument, NULL, 'b' },
	{ "threads",          required_argument, NULL, 'j' },
	{ "kmer",             required_argument, NULL, 'k' },
	{ "chastity",         no_argument, &opt::chastityFilter, 1 },
	{ "no-chastity",      no_argument, &opt::chastityFilter, 0 },
	{ "trim-masked",      no_argument, &opt::trimMasked, 1 },
	{ "no-trim-masked",   no_argument, &opt::trimMasked, 0 },
	{ "trim-quality",     required_argument, NULL, 'q' },
	{ "standard-quality", no_argument, &opt::qualityOffset, 33 },
	{ "illumina-quality", no_argument, &opt::qualityOffset, 64 },
	{ "verbose",          no_argument, NULL, 'v' },
	{ "help",             no_argument, NULL, OPT_HELP },
	{ "version",          no_argument, NULL, OPT_VERSION },
	{ NULL, 0, NULL, 0 }
};

template <typename It>
void loadBloomFilter(DBGBloom& g, It first, It last)
{
	if (opt::verbose > 0)
		cerr << "Constructing the Bloom filter\n";
	for (It it = first; it != last; ++it) {
		std::string path = *it;
		if (opt::verbose > 0)
			cerr << "Reading `" << path << "'...\n";
		g.open(path, opt::verbose >= 2);
	}
	if (opt::verbose > 0)
		cerr << "Loaded " << num_vertices(g) << " k-mer\n";
}

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
		  case 'j':
			arg >> opt::threads; break;
		  case 'k':
			arg >> opt::k; break;
		  case 'q':
			arg >> opt::qualityThreshold; break;
		  case 'v':
			opt::verbose++; break;
		  case OPT_HELP:
			cerr << USAGE_MESSAGE;
			exit(EXIT_SUCCESS);
		  case OPT_VERSION:
			cerr << VERSION_MESSAGE;
			exit(EXIT_SUCCESS);
		}
		if (optarg != NULL && (!arg.eof() || arg.fail())) {
			cerr << PROGRAM ": invalid option: `-"
				<< (char)c << optarg << "'\n";
			exit(EXIT_FAILURE);
		}
	}

	if (opt::k == 0) {
		cerr << PROGRAM ": missing mandatory option `-k'\n";
		die = true;
	}

	if (argc - optind < 2) {
		cerr << PROGRAM ": missing arguments\n";
		die = true;
	}

	if (die) {
		cerr << "Try `" << PROGRAM
			<< " --help' for more information.\n";
		exit(EXIT_FAILURE);
	}

	assert_good(cout, "stdout");

	return 0;
}
