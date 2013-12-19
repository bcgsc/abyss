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
"Usage 1: " PROGRAM " build [GLOBAL_OPTS] [COMMAND_OPTS] <OUTPUT_BLOOM_FILE> <READS_FILE_1> [<READS_FILE_2>]...\n"
"Usage 2: " PROGRAM " union [GLOBAL_OPTS] [COMMAND_OPTS] <OUTPUT_BLOOM_FILE> <BLOOM_FILE_1> [<BLOOM_FILE_2>]...\n"
"Build and manipulate bloom filter files.\n"
"\n"
" Global options:\n"
"\n"
"  -k, --kmer=N               the size of a k-mer [required]\n"
"  -v, --verbose              display verbose output\n"
"      --help                 display this help and exit\n"
"      --version              output version information and exit\n"
"\n"
" Options for `" PROGRAM " build':\n"
"\n"
"  -b, --bloom-size=N         size of bloom filter [500M]\n"
"  -l, --levels=N             build a counting bloom filter with N levels\n"
"                             and output the last level\n"
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
"\n"
" Options for `" PROGRAM " union': (none)\n"
"\n"
"Report bugs to <" PACKAGE_BUGREPORT ">.\n";

namespace opt {
	/** The size of the bloom filter in bytes. */
	size_t bloomSize = 500 * 1024 * 1024;

	/** The size of a k-mer. */
	unsigned k;

	/** Number of levels for counting bloom filter. */
	unsigned levels = 1;
}

static const char shortopts[] = "b:k:l:q:v";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
	{ "bloom-size",       required_argument, NULL, 'b' },
	{ "kmer",             required_argument, NULL, 'k' },
	{ "levels",           required_argument, NULL, 'l' },
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

void parseGlobalOpts(int argc, char** argv)
{
	bool done = false;
	int optindPrev = optind;

	for (int c; (c = getopt_long(argc, argv,
					shortopts, longopts, NULL)) != -1;) {

		istringstream arg(optarg != NULL ? optarg : "");
		switch (c) {
		  case '?':
			dieWithUsageError();
		  case 'k':
			arg >> opt::k; break;
		  case 'v':
			opt::verbose++; break;
		  case OPT_HELP:
			cerr << USAGE_MESSAGE;
			exit(EXIT_SUCCESS);
		  case OPT_VERSION:
			cerr << VERSION_MESSAGE;
			exit(EXIT_SUCCESS);
		  default:
			// end of global opts
			optind = optindPrev;
			done = true;
			break;
		}

		if (done)
			break;

		if (optarg != NULL && (!arg.eof() || arg.fail())) {
			cerr << PROGRAM ": invalid option: `-"
				<< (char)c << optarg << "'\n";
			exit(EXIT_FAILURE);
		}

		optindPrev = optind;
	}

	if (opt::k == 0) {
		cerr << PROGRAM ": missing mandatory option `-k'\n";
		dieWithUsageError();
	}

	Kmer::setLength(opt::k);
}

int build(int argc, char** argv)
{
	parseGlobalOpts(argc, argv);

	for (int c; (c = getopt_long(argc, argv,
					shortopts, longopts, NULL)) != -1;) {
		istringstream arg(optarg != NULL ? optarg : "");
		switch (c) {
		  case '?':
			dieWithUsageError();
		  case 'b':
			opt::bloomSize = SIToBytes(arg); break;
		  case 'l':
			arg >> opt::levels; break;
		  case 'q':
			arg >> opt::qualityThreshold; break;
		}
		if (optarg != NULL && (!arg.eof() || arg.fail())) {
			cerr << PROGRAM ": invalid option: `-"
				<< (char)c << optarg << "'\n";
			exit(EXIT_FAILURE);
		}
	}

	if (opt::levels > 2)
	{
		cerr << PROGRAM ": --levels > 2 is not currently supported\n";
		dieWithUsageError();
	}

	if (argc - optind < 2) {
		cerr << PROGRAM ": missing arguments\n";
		dieWithUsageError();
	}

	char* outputPath = argv[optind];
	optind++;

	if (opt::levels == 1) {

		// Specify bloom filter size in bits.
		BloomFilter bloom(opt::bloomSize * 8);

		for (int i = optind; i < argc; i++)
			bloom.loadFile(opt::k, argv[i], opt::verbose);

		if (opt::verbose)
			cerr << "Writing bloom filter to `"
				<< outputPath << "'...\n";

		ofstream outputFile(outputPath, ios_base::out | ios_base::binary);
		assert_good(outputFile, outputPath);
		outputFile << bloom;
		outputFile.flush();
		assert_good(outputFile, outputPath);
		outputFile.close();


	} else {

		// Specify bloom filter size in bits. Divide by two
		// because counting bloom filter requires twice as
		// much space.
		CountingBloomFilter bloom(opt::bloomSize * 8 / 2);

		for (int i = optind; i < argc; i++)
			bloom.loadFile(opt::k, argv[i], opt::verbose);

		if (opt::verbose)
			cerr << "Writing bloom filter to `"
				<< outputPath << "'...\n";

		ofstream outputFile(outputPath, ios_base::out | ios_base::binary);
		assert_good(outputFile, outputPath);
		outputFile << bloom.getBloomFilter(1);
		outputFile.flush();
		assert_good(outputFile, outputPath);
		outputFile.close();

	}

	return 0;
}

int union_(int argc, char** argv)
{
	parseGlobalOpts(argc, argv);

	if (argc - optind < 3) {
		cerr << PROGRAM ": missing arguments\n";
		dieWithUsageError();
	}

	const char* outputPath = argv[optind];
	optind++;

	BloomFilter unionBloom;

	for (int i = optind; i < argc; i++) {
		const char* path = argv[i];
		if (opt::verbose)
			std::cerr << "Loading bloom filter from `"
				<< path << "'...\n";
		ifstream input(path, ios_base::in | ios_base::binary);
		assert_good(input, path);
		// The second arg is the "load as union" flag.
		// Loading the first input BF with the flag as false
		// sets the required size for remaining input bloom filters.
		unionBloom.read(input, i > optind);
		assert_good(input, path);
		input.close();
	}

	if (opt::verbose)
		std::cerr << "Writing union of bloom filters to `"
			<< outputPath << "'...\n";
	ofstream output(outputPath, ios_base::out | ios_base::binary);
	assert_good(output, outputPath);
	output << unionBloom;
	output.flush();
	assert_good(output, outputPath);
	output.close();

	return 0;
}

int main(int argc, char** argv)
{
	if (argc < 2)
		dieWithUsageError();

	string command(argv[1]);
	optind++;

	if (command == "build")
		return build(argc, argv);
	else if (command == "union")
		return union_(argc, argv);

	dieWithUsageError();
}
