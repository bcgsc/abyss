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

#if _OPENMP
# include <omp.h>
# include "Bloom/ConcurrentBloomFilter.h"
#endif

using namespace std;

#define PROGRAM "abyss-bloom"

static const char VERSION_MESSAGE[] =
PROGRAM " (" PACKAGE_NAME ") " VERSION "\n"
"Written by Shaun Jackman, Hamid Mohamadi, Anthony Raymond and\n"
"Ben Vandervalk.\n"
"\n"
"Copyright 2013 Canada's Michael Smith Genome Science Centre\n";

static const char USAGE_MESSAGE[] =
"Usage 1: " PROGRAM " build [GLOBAL_OPTS] [COMMAND_OPTS] <OUTPUT_BLOOM_FILE> <READS_FILE_1> [READS_FILE_2]...\n"
"Usage 2: " PROGRAM " union [GLOBAL_OPTS] [COMMAND_OPTS] <OUTPUT_BLOOM_FILE> <BLOOM_FILE_1> <BLOOM_FILE_2> [BLOOM_FILE_3]...\n"
"Usage 3: " PROGRAM " intersect [GLOBAL_OPTS] [COMMAND_OPTS] <OUTPUT_BLOOM_FILE> <BLOOM_FILE_1> <BLOOM_FILE_2> [BLOOM_FILE_3]...\n"
"Usage 4: " PROGRAM " info [GLOBAL_OPTS] [COMMAND_OPTS] <BLOOM_FILE>\n"
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
"  -j, --threads=N            use N parallel threads [1]\n"
"  -l, --levels=N             build a counting bloom filter with N levels\n"
"                             and output the last level\n"
"  -L, --init-level='N=FILE'  initialize level N of counting bloom filter\n"
"                             from FILE\n"
"      --chastity             discard unchaste reads [default]\n"
"      --no-chastity          do not discard unchaste reads\n"
"      --trim-masked          trim masked bases from the ends of reads\n"
"      --no-trim-masked       do not trim masked bases from the ends\n"
"                             of reads [default]\n"
"  -n, --num-locks=N          number of write locks on bloom filter\n"
"  -q, --trim-quality=N       trim bases from the ends of reads whose\n"
"                             quality is less than the threshold\n"
"      --standard-quality     zero quality is `!' (33)\n"
"                             default for FASTQ and SAM files\n"
"      --illumina-quality     zero quality is `@' (64)\n"
"                             default for qseq and export files\n"
"  -w, --window M/N           build a bloom filter for subwindow M of N\n"
"\n"
" Options for `" PROGRAM " union': (none)\n"
" Options for `" PROGRAM " intersect': (none)\n"
" Options for `" PROGRAM " info': (none)\n"
"\n"
"Report bugs to <" PACKAGE_BUGREPORT ">.\n";

namespace opt {
	/** The size of the bloom filter in bytes. */
	size_t bloomSize = 500 * 1024 * 1024;

	/** The number of parallel threads. */
	unsigned threads = 1;

	/** The size of a k-mer. */
	unsigned k;

	/** Number of levels for counting bloom filter. */
	unsigned levels = 1;

	/**
	 * Files used to initialize levels of counting
	 * bloom filter (-L option).
	 */
	vector< vector<string> > levelInitPaths;

	/**
	 * Num of locked windows to use, when invoking with
	 * the -j option.
	 */
	size_t numLocks = 1000;

	/** Index of bloom filter window.
	  ("M" for -w option) */
	unsigned windowIndex = 0;

	/** Number of windows in complete bloom filter.
	  ("N" for -w option) */
	unsigned windows = 0;
}

static const char shortopts[] = "b:j:k:l:L:n:q:vw:";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
	{ "bloom-size",       required_argument, NULL, 'b' },
	{ "threads",          required_argument, NULL, 'j' },
	{ "kmer",             required_argument, NULL, 'k' },
	{ "levels",           required_argument, NULL, 'l' },
	{ "init-level",       required_argument, NULL, 'L' },
	{ "chastity",         no_argument, &opt::chastityFilter, 1 },
	{ "no-chastity",      no_argument, &opt::chastityFilter, 0 },
	{ "trim-masked",      no_argument, &opt::trimMasked, 1 },
	{ "no-trim-masked",   no_argument, &opt::trimMasked, 0 },
	{ "num-locks",        required_argument, NULL, 'n' },
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
			cout << USAGE_MESSAGE;
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
			Bloom::LoadType loadType = (j > 0) ?
				Bloom::LOAD_UNION : Bloom::LOAD_OVERWRITE;
			bf.getBloomFilter(i).read(*in, loadType);
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

	if (opt::verbose)
		cerr << "Successfully loaded bloom filter.\n";
}

template <typename BF>
void writeBloom(BF& bf, string& outputPath)
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

template <typename BF>
void printBloomStats(ostream& os, const BF& bloom)
{
	os << "Bloom size (bits): " << bloom.size() << "\n"
		<< "Bloom popcount (bits): " << bloom.popcount() << "\n"
		<< "Bloom filter FPR: " << setprecision(3)
			<< 100 * bloom.FPR() << "%\n";
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
		  case 'j':
			arg >> opt::threads; break;
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
		  case 'n':
			arg >> opt::numLocks; break;
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

#if _OPENMP
	if (opt::threads > 0)
		omp_set_num_threads(opt::threads);
#endif

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
#ifdef _OPENMP
			ConcurrentBloomFilter<BloomFilter>
				cbf(bloom, opt::numLocks);
			loadFilters(cbf, argc, argv);
#else
			loadFilters(bloom, argc, argv);
#endif
			printBloomStats(cerr, bloom);
			writeBloom(bloom, outputPath);
		}
		else {
			CascadingBloomFilter cascadingBloom(bits);
			initBloomFilterLevels(cascadingBloom);
#ifdef _OPENMP
			ConcurrentBloomFilter<CascadingBloomFilter>
				cbf(cascadingBloom, opt::numLocks);
			loadFilters(cbf, argc, argv);
#else
			loadFilters(cascadingBloom, argc, argv);
#endif
			printBloomStats(cerr, cascadingBloom);
			writeBloom(cascadingBloom, outputPath);
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
			printBloomStats(cerr, bloom);
			writeBloom(bloom, outputPath);
		}
		else {
			CascadingBloomFilterWindow cascadingBloom(bits, startBitPos, endBitPos);
			initBloomFilterLevels(cascadingBloom);
			loadFilters(cascadingBloom, argc, argv);
			printBloomStats(cerr, cascadingBloom);
			writeBloom(cascadingBloom, outputPath);
		}
	}

	return 0;
}

int combine(int argc, char** argv, Bloom::LoadType loadType)
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
		Bloom::LoadType loadOp = (i > optind) ?
				loadType : Bloom::LOAD_OVERWRITE;
		bloom.read(*in, loadOp);
		assert_good(*in, path);
		closeInputStream(in, path);
	}

	if (opt::verbose) {
		cerr << "Successfully loaded bloom filter.\n";
		printBloomStats(cerr, bloom);
		switch(loadType) {
			case Bloom::LOAD_UNION:
				std::cerr << "Writing union of bloom filters to `"
					<< outputPath << "'...\n";
				break;
			case Bloom::LOAD_INTERSECT:
				std::cerr << "Writing intersection of bloom filters to `"
					<< outputPath << "'...\n";
				break;
			default:
				std::cerr << "error: expected LOAD_UNION or LOAD_INTERSECT\n";
				assert(false);
				break;
		}
	}

	ostream* out = openOutputStream(outputPath);

	assert_good(*out, outputPath);
	*out << bloom;
	out->flush();
	assert_good(*out, outputPath);

	closeOutputStream(out, outputPath);

	return 0;
}

int info(int argc, char** argv)
{
	parseGlobalOpts(argc, argv);

	if (argc - optind < 1) {
		cerr << PROGRAM ": missing arguments\n";
		dieWithUsageError();
	}

	BloomFilter bloom;
	string path = argv[optind];

	if (opt::verbose)
		std::cerr << "Loading bloom filter from `"
			<< path << "'...\n";

	istream* in = openInputStream(path);
	assert_good(*in, path);
	*in >> bloom;

	printBloomStats(cerr, bloom);

	closeInputStream(in, path);

	return 0;
}

int main(int argc, char** argv)
{
	if (argc < 2)
		dieWithUsageError();

	string command(argv[1]);
	optind++;

	if (command == "--help" || command == "-h") {
		cout << USAGE_MESSAGE;
		return EXIT_SUCCESS;
	}
	else if (command == "build") {
		return build(argc, argv);
	}
	else if (command == "union") {
		return combine(argc, argv, Bloom::LOAD_UNION);
	}
	else if (command == "intersect") {
		return combine(argc, argv, Bloom::LOAD_INTERSECT);
	}
	else if (command == "info") {
		return info(argc, argv);
	}

	dieWithUsageError();
}
