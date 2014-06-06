/**
 * Build and manipulate bloom filter files.
 */

#include "config.h"
#include "Common/Options.h"
#include "Common/Kmer.h"
#include "DataLayer/Options.h"
#include "Common/StringUtil.h"
#include "Bloom/Bloom.h"
#include "Bloom/BloomFilter.h"
#include "Bloom/CascadingBloomFilter.h"
#include "Bloom/BloomFilterWindow.h"
#include "Bloom/CascadingBloomFilterWindow.h"

#include <cstdlib>
#include <getopt.h>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

#define PROGRAM "abyss-bloom"

static const char VERSION_MESSAGE[] =
PROGRAM " (" PACKAGE_NAME ") " VERSION "\n"
"Written by Shaun Jackman, Hamid Mohamadi, Anthony Raymond, \n"
"Ben Vandervalk, and Justin Chu.\n"
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
"  -L, --init-level='N=FILE'  initialize level N of counting bloom filter\n"
"                             from FILE\n"
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
"  -w, --window M/N           build a bloom filter for subwindow M of N\n"
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

	/**
	 * Files used to initialize levels of counting
	 * bloom filter (-L option).
	 */
	vector< vector<string> > levelInitPaths;

	/** Index of bloom filter window.
	  ("M" for -w option) */
	unsigned windowIndex = 0;

	/** Number of windows in complete bloom filter.
	  ("N" for -w option) */
	unsigned windows = 0;
}

static const char shortopts[] = "b:k:l:L:q:vw:";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
	{ "bloom-size",       required_argument, NULL, 'b' },
	{ "kmer",             required_argument, NULL, 'k' },
	{ "levels",           required_argument, NULL, 'l' },
	{ "init-level",       required_argument, NULL, 'L' },
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
	{ "window",           required_argument, NULL, 'w' },
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

static inline istream* openInputStream(const string& path)
{
	if (path == "-")
		return &cin;
	return new ifstream(path.c_str());
}

static inline ostream* openOutputStream(const string& path)
{
	if (path == "-")
		return &cout;
	return new ofstream(path.c_str());
}

static inline void closeInputStream(istream* in, const string& path)
{
	if (path == "-")
		return;
	ifstream* ifs = static_cast<ifstream*>(in);
	ifs->close();
	delete ifs;
}

static inline void closeOutputStream(ostream* out, const string& path)
{
	if (path == "-")
		return;
	ofstream* ofs = static_cast<ofstream*>(out);
	ofs->close();
	delete ofs;
}

template <typename CBF>
void initBloomFilterLevels(CBF& bf)
{
	assert(opt::levels >= 2);
	assert(opt::levelInitPaths.size() <= opt::levels);

	for (unsigned i = 0; i < opt::levelInitPaths.size(); i++) {
		vector<string>& paths = opt::levelInitPaths.at(i);
		for (unsigned j = 0; j < paths.size(); j++) {
			string path = paths.at(j);
			cerr << "Loading `" << path << "' into level "
				<< i + 1 << " of counting bloom filter...\n";
			istream* in = openInputStream(path);
			assert(*in);
			bf.getBloomFilter(i).read(*in, j > 0);
			assert(*in);
			closeInputStream(in, path);
		}
	}
}

template <typename BF>
void loadFilters(BF& bf, int argc, char** argv)
{
	for (int i = optind; i < argc; i++)
		Bloom::loadFile(bf, opt::k, argv[i], opt::verbose);

	if (opt::verbose) {
		cerr << "Successfully loaded bloom filter.\n";
		cerr << "Bloom filter FPR: " << setprecision(3)
			<< 100 * bf.FPR() << "%\n";
	}
}

template <typename BF>
void buildAndOutput(BF& bf, string& outputPath)
{
	if (opt::verbose) {
		cerr << "Writing bloom filter to `"
			<< outputPath << "'...\n";
	}

	ostream* out = openOutputStream(outputPath);

	assert_good(*out, outputPath);
	*out << bf;
	out->flush();
	assert_good(*out, outputPath);

	closeOutputStream(out, outputPath);
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
		  case 'L':
			{
				unsigned level;
				arg >> level >> expect("=");
				if (arg.fail() || arg.eof())
					break;
				string path;
				arg >> path;
				if (arg.fail())
					break;
				if (level > opt::levelInitPaths.size())
					opt::levelInitPaths.resize(level);
				opt::levelInitPaths[level-1].push_back(path);
				break;
			}
		  case 'q':
			arg >> opt::qualityThreshold; break;
		  case 'w':
			arg >> opt::windowIndex;
			arg >> expect("/");
			arg >> opt::windows;
			break;
		}
		if (optarg != NULL && (!arg.eof() || arg.fail())) {
			cerr << PROGRAM ": invalid option: `-"
				<< (char)c << optarg << "'\n";
			exit(EXIT_FAILURE);
		}
	}

	if (opt::levels > 2)
	{
		cerr << PROGRAM ": -l > 2 is not currently supported\n";
		dieWithUsageError();
	}

	if (!opt::levelInitPaths.empty() && opt::levels < 2)
	{
		cerr << PROGRAM ": -L can only be used with counting bloom "
			"filters (-l >= 2)\n";
		dieWithUsageError();
	}

	if (opt::levelInitPaths.size() > opt::levels) {
		cerr << PROGRAM ": level arg to -L is greater than number"
			" of bloom filter levels (-l)\n";
		dieWithUsageError();
	}

	// bloom filter size in bits
	size_t bits = opt::bloomSize * 8;

	if (bits % opt::levels != 0) {
		cerr << PROGRAM ": bloom filter size (-b) must be evenly divisible "
			<< "by number of bloom filter levels (-l)\n";
		dieWithUsageError();
	}

	if (opt::windows != 0 && bits / opt::levels % opt::windows != 0) {
		cerr << PROGRAM ": (b / l) % w == 0 must be true, where "
			<< "b is bloom filter size (-b), "
			<< "l is number of levels (-l), and "
			<< "w is number of windows (-w)\n";
		dieWithUsageError();
	}

	if (argc - optind < 2) {
		cerr << PROGRAM ": missing arguments\n";
		dieWithUsageError();
	}

	// if we are building a counting bloom filter, reduce
	// the size of each level so that the overall bloom filter
	// fits within the memory limit (specified by -b)
	bits /= opt::levels;

	string outputPath(argv[optind]);
	optind++;

	if (opt::windows == 0) {

		if (opt::levels == 1) {
			BloomFilter bloom(bits);
			loadFilters(bloom, argc, argv);
			buildAndOutput(bloom, outputPath);
		}
		else {
			CascadingBloomFilter cascadeBloom(bits);
			initBloomFilterLevels(cascadeBloom);
			loadFilters(cascadeBloom, argc, argv);
			buildAndOutput(cascadeBloom, outputPath);
		}

	} else {

		size_t bitsPerWindow = bits / opt::windows;
		size_t startBitPos = (opt::windowIndex - 1) * bitsPerWindow;
		size_t endBitPos;

		if (opt::windowIndex < opt::windows)
			endBitPos = opt::windowIndex * bitsPerWindow - 1;
		else
			endBitPos = bits - 1;

		if (opt::levels == 1) {
			BloomFilterWindow bloom(bits, startBitPos, endBitPos);
			loadFilters(bloom, argc, argv);
			buildAndOutput(bloom, outputPath);
		} else {
			CascadingBloomFilterWindow cascadeBloom(bits, startBitPos,
					endBitPos);
			initBloomFilterLevels(cascadeBloom);
			loadFilters(cascadeBloom, argc, argv);
			buildAndOutput(cascadeBloom, outputPath);
		}

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

	string outputPath(argv[optind]);
	optind++;

	BloomFilter bloom;

	for (int i = optind; i < argc; i++) {
		string path(argv[i]);
		if (opt::verbose)
			std::cerr << "Loading bloom filter from `"
				<< path << "'...\n";
		istream* in = openInputStream(path);
		assert_good(*in, path);
		bloom.read(*in, i > optind);
		assert_good(*in, path);
		closeInputStream(in, path);
	}

	if (opt::verbose) {
		cerr << "Successfully loaded bloom filter.\n";
		cerr << "Bloom filter FPR: " << setprecision(3)
			<< 100 * bloom.FPR() << "%\n";
	}
	buildAndOutput(bloom, outputPath);

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
