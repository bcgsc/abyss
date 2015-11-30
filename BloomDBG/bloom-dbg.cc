#include "config.h"

#include "Common/Kmer.h"

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
	"Anthony Raymond, and Justin Chu.\n"
	"\n"
	"Copyright 2015 Canada's Michael Smith Genome Science Centre\n";

static const char USAGE_MESSAGE[] =
"Usage: " PROGRAM " -k <kmer_size> [options]... <FASTQ> [FASTQ]... > assembly.fasta\n"
"Perform a de Bruijn graph assembly of the given FASTQ files.\n"
"\n"
" Options:\n"
"\n"
"      --help                 display this help and exit\n"
"  -j, --threads=N            use N parallel threads [1]\n"
"  -k, --kmer=N               the size of a k-mer\n"
"      --version              output version information and exit\n"
"\n"
"Report bugs to <" PACKAGE_BUGREPORT ">.\n";

namespace opt {

	/** The number of parallel threads. */
	static unsigned threads = 1;

	/** The size of a k-mer. */
	unsigned k;

}

static const char shortopts[] = "j:k:";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
	{ "help",             no_argument, NULL, OPT_HELP },
	{ "threads",          required_argument, NULL, 'j' },
	{ "kmer",             required_argument, NULL, 'k' },
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
		  case 'j':
			arg >> opt::threads; break;
		  case 'k':
			arg >> opt::k; break;
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

	Kmer::setLength(opt::k);

	return EXIT_SUCCESS;
}
