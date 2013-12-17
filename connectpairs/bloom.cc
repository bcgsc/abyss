/**
 * Build and manipulate bloom filter files.
 */

#include "config.h"
#include "Common/Options.h"
#include "Common/Kmer.h"
#include "DataLayer/Options.h"
#include "Common/StringUtil.h"
#include "connectpairs/CountingBloomFilter.h"

#include <cstdlib>
#include <getopt.h>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

#define PROGRAM "abyss-bloom"

static const char VERSION_MESSAGE[] =
PROGRAM " (" PACKAGE_NAME ") " VERSION "\n"
"Written by Shaun Jackman, Hamid Mohamadi, Anthony Raymond and\n"
"Ben Vandervalk.\n"
"\n"
"Copyright 2013 Canada's Michael Smith Genome Science Centre\n";

static const char USAGE_MESSAGE[] =
"Usage: " PROGRAM "build [OPTIONS] <OUTPUT_BLOOM_FILE> <READS_FILE_1> [<READS_FILE_2>]...\n"
"Build and manipulate bloom filter files.\n"
"\n"
" Options:\n"
"\n"
"  -b, --bloom-size=N         size of bloom filter [500M]\n"
"  -k, --kmer=N               the size of a k-mer\n"
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

namespace opt {
	/** The size of the bloom filter in bytes. */
	size_t bloomSize = 500 * 1024 * 1024;

	/** The size of a k-mer. */
	unsigned k;
}

static const char shortopts[] = "b:k:q:v";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
	{ "bloom-size",       required_argument, NULL, 'b' },
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

void dieWithUsageError()
{
	cerr << "Try `" << PROGRAM
		<< " --help' for more information.\n";
	exit(EXIT_FAILURE);
}

int main(int argc, char** argv)
{
	if (argc < 2)
		dieWithUsageError();

	string command(argv[1]);

	if (command != "build")
		dieWithUsageError();

	optind++;

	for (int c; (c = getopt_long(argc, argv,
					shortopts, longopts, NULL)) != -1;) {
		istringstream arg(optarg != NULL ? optarg : "");
		switch (c) {
		  case '?':
			dieWithUsageError();
		  case 'b':
			opt::bloomSize = SIToBytes(arg); break;
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
		dieWithUsageError();
	}

	if (argc - optind < 2) {
		cerr << PROGRAM ": missing arguments\n";
		dieWithUsageError();
	}

	Kmer::setLength(opt::k);

	// Specify bloom filter size in bits. Divide by two
	// because counting bloom filter requires twice as
	// much space.
	size_t bits = opt::bloomSize * 8 / 2;
	CountingBloomFilter bloom(bits);

	for (int i = optind + 1; i < argc; i++)
		bloom.loadFile(opt::k, argv[i], opt::verbose);

	char* outputPath = argv[optind];
	if (opt::verbose)
		cerr << "Writing bloom filter to `"
			<< outputPath << "'...\n";

	ofstream outputFile(outputPath, ios_base::out | ios_base::binary);
	assert_good(outputFile, outputPath);
	outputFile << bloom.getBloomFilter(1);
	outputFile.flush();
	assert_good(outputFile, outputPath);
	outputFile.close();

	return 0;
}
