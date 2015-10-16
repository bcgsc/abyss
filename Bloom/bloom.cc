/**
 * Build and manipulate bloom filter files.
 */

#include "config.h"
#include "Common/Options.h"
#include "Common/Kmer.h"
#include "Common/BitUtil.h"
#include "Common/KmerIterator.h"
#include "Graph/Path.h"
#include "Graph/ExtendPath.h"
#include "Konnector/DBGBloom.h"
#include "DataLayer/Options.h"
#include "DataLayer/FastaReader.h"
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
#include <cmath>

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
"Usage 5: " PROGRAM " compare [GLOBAL_OPTS] [COMMAND_OPTS] <BLOOM_FILE_1> <BLOOM_FILE_2>\n"
"Usage 6: " PROGRAM " kmers [GLOBAL_OPTS] [COMMAND_OPTS] <BLOOM_FILE> <READS_FILE>\n"
"Usage 7: " PROGRAM " trim [GLOBAL_OPTS] [COMMAND_OPTS] <BLOOM_FILE> <READS_FILE> [READS_FILE_2]... > trimmed.fq\n"
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
"  -B, --buffer-size=N        size of I/O buffer for each thread, in bytes [100000]\n"
"  -j, --threads=N            use N parallel threads [1]\n"
"  -l, --levels=N             build a cascading bloom filter with N levels\n"
"                             and output the last level\n"
"  -L, --init-level='N=FILE'  initialize level N of cascading bloom filter\n"
"                             from FILE\n"
"      --chastity             discard unchaste reads [default]\n"
"      --no-chastity          do not discard unchaste reads\n"
"      --trim-masked          trim masked bases from the ends of reads\n"
"      --no-trim-masked       do not trim masked bases from the ends\n"
"                             of reads [default]\n"
"  -n, --num-locks=N          number of write locks on bloom filter [1000]\n"
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
" Options for `" PROGRAM " compare':\n"
"\n"
"  -m, --method=`String'      choose distance calculation method \n"
"                             [`jaccard'(default), `forbes', `czekanowski']\n"
"\n"
" Options for `" PROGRAM " kmers':\n"
"\n"
"  -r, --inverse              get k-mers that are *NOT* in the bloom filter\n"
"  --bed                      output k-mers in BED format\n"
"  --fasta                    output k-mers in FASTA format [default]\n"
"  --raw                      output k-mers in raw format (one per line)\n"
"\n"
" Options for `" PROGRAM " trim': (none)\n"
"\n"
"Report bugs to <" PACKAGE_BUGREPORT ">.\n";;

enum OutputFormat { BED, FASTA, RAW };

namespace opt {

	/** The size of the bloom filter in bytes. */
	size_t bloomSize = 500 * 1024 * 1024;

	/** The size of the I/O buffer of each thread, in bytes  */
	size_t bufferSize = 100000;

	/** The number of parallel threads. */
	unsigned threads = 1;

	/** The size of a k-mer. */
	unsigned k;

	/** Number of levels for cascading bloom filter. */
	unsigned levels = 1;

	/**
	 * Files used to initialize levels of cascading
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

	/* Method for similarity or distance calculation.
	 -m option
	 */
	string method("jaccard");

	/* Inverse option to retrieve kmers which are not
	 in the filter
	 */
	bool inverse = false;

	OutputFormat format = FASTA;
}

static const char shortopts[] = "b:B:j:k:l:L:m:n:q:rvw:";

enum { OPT_HELP = 1, OPT_VERSION, OPT_BED, OPT_FASTA, OPT_RAW };

static const struct option longopts[] = {
	{ "bloom-size", required_argument, NULL, 'b' },
	{ "buffer-size", required_argument, NULL, 'B' },
	{ "threads", required_argument, NULL, 'j' },
	{ "kmer", required_argument, NULL, 'k' },
	{ "levels", required_argument, NULL, 'l' },
	{ "init-level", required_argument, NULL, 'L' },
	{ "chastity", no_argument, &opt::chastityFilter, 1 },
	{ "no-chastity", no_argument, &opt::chastityFilter, 0 },
	{ "trim-masked", no_argument, &opt::trimMasked, 1 },
	{ "no-trim-masked", no_argument, &opt::trimMasked, 0 },
	{ "num-locks", required_argument, NULL, 'n' },
	{ "trim-quality", required_argument, NULL, 'q' },
	{ "standard-quality", no_argument, &opt::qualityOffset, 33 },
	{ "illumina-quality", no_argument, &opt::qualityOffset, 64 },
	{ "verbose", no_argument, NULL, 'v' },
	{ "help", no_argument, NULL, OPT_HELP },
	{ "version", no_argument, NULL, OPT_VERSION },
	{ "window", required_argument, NULL, 'w' },
	{ "method", required_argument, NULL, 'm' },
	{ "inverse", required_argument, NULL, 'r' },
	{ "bed", no_argument, NULL, OPT_BED },
	{ "fasta", no_argument, NULL, OPT_FASTA },
	{ "raw", no_argument, NULL, OPT_RAW },
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
			cout << VERSION_MESSAGE;
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
				<< i + 1 << " of cascading bloom filter...\n";
			istream* in = openInputStream(path);
			assert(*in);
			BitwiseOp readOp = (j > 0) ? BITWISE_OR : BITWISE_OVERWRITE;
			bf.getBloomFilter(i).read(*in, readOp);
			assert(*in);
			closeInputStream(in, path);
		}
	}
}

template <typename BF>
void loadFilters(BF& bf, int argc, char** argv)
{
	for (int i = optind; i < argc; i++)
		Bloom::loadFile(bf, opt::k, argv[i], opt::verbose, opt::bufferSize);

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

template <typename BF>
void printCascadingBloomStats(ostream& os, BF& bloom)
{
	for (unsigned i = 0; i < opt::levels; i++) {
		os << "Stats for Bloom filter level " << i+1 << ":\n"
			<< "\tBloom size (bits): "
			<< bloom.getBloomFilter(i).size() << "\n"
			<< "\tBloom popcount (bits): "
			<< bloom.getBloomFilter(i).popcount() << "\n"
			<< "\tBloom filter FPR: " << setprecision(3)
			<< 100 * bloom.getBloomFilter(i).FPR() << "%\n";
	}
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
		  case 'B':
			arg >> opt::bufferSize; break;
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

	if (!opt::levelInitPaths.empty() && opt::levels < 2)
	{
		cerr << PROGRAM ": -L can only be used with cascading bloom "
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

	// if we are building a cascading bloom filter, reduce
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
			CascadingBloomFilter cascadingBloom(bits, opt::levels);
			initBloomFilterLevels(cascadingBloom);
#ifdef _OPENMP
			ConcurrentBloomFilter<CascadingBloomFilter>
				cbf(cascadingBloom, opt::numLocks);
			loadFilters(cbf, argc, argv);
#else
			loadFilters(cascadingBloom, argc, argv);
#endif
			printCascadingBloomStats(cerr, cascadingBloom);
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
			CascadingBloomFilterWindow cascadingBloom(
				bits, startBitPos, endBitPos, opt::levels);
			initBloomFilterLevels(cascadingBloom);
			loadFilters(cascadingBloom, argc, argv);
			printCascadingBloomStats(cerr, cascadingBloom);
			writeBloom(cascadingBloom, outputPath);
		}
	}

	return 0;
}

int combine(int argc, char** argv, BitwiseOp readOp)
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
		BitwiseOp op = (i > optind) ? readOp : BITWISE_OVERWRITE;
		bloom.read(*in, op);
		assert_good(*in, path);
		closeInputStream(in, path);
	}

	if (opt::verbose) {
		cerr << "Successfully loaded bloom filter.\n";
		printBloomStats(cerr, bloom);
		switch(readOp) {
			case BITWISE_OR:
				std::cerr << "Writing union of bloom filters to `"
					<< outputPath << "'...\n";
				break;
			case BITWISE_AND:
				std::cerr << "Writing intersection of bloom filters to `"
					<< outputPath << "'...\n";
				break;
			default:
				std::cerr << "error: expected BITWISE_OR or BITWISE_AND\n";
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

int compare(int argc, char ** argv){
	parseGlobalOpts(argc, argv);
	// Arg parser to get `m' option in case set
	for (int c; (c = getopt_long(argc, argv,
								 shortopts, longopts, NULL)) != -1;) {
		istringstream arg(optarg != NULL ? optarg : "");
		switch (c) {
			case '?':
				cerr << PROGRAM ": unrecognized option: `-" << optopt
					<< "'" << endl;
				dieWithUsageError();
			case 'm':
				arg >> opt::method; break;
				break;
		}
		if (optarg != NULL && (!arg.eof() || arg.fail())) {
			cerr << PROGRAM ": invalid option: `-"
			<< (char)c << optarg << "'\n";
			exit(EXIT_FAILURE);
		}
		if (opt::method != "jaccard" && opt::method != "czekanowski" && opt::method != "forbes")
			std::cerr << "Invalid method: " << opt::method << std::endl;
	}


	// Set method strin
	string method(opt::method);
	if (opt::verbose)
	std::cerr << "Computing distance for 2"
			  << " samples...\n";
	// Get both paths and open istreams
	BloomFilter bloomA;
	string pathA(argv[optind]);
	BloomFilter bloomB;
	string pathB(argv[optind+1]);
	if (opt::verbose)
	  std::cerr << "Loading bloom filters from "
		<< pathA << " and " << pathB << "...\n";
	istream* inA = openInputStream(pathA);
	istream* inB = openInputStream(pathB);
	// Assert state of streams
	assert_good(*inA, pathA);
	assert_good(*inB, pathB);
	// Not sure this conversion is needed, check docs
	std::istream & tA = *inA;
	std::istream & tB = *inB;
	// Need to read header for bit start and end info
	Bloom::FileHeader headerA = Bloom::readHeader(tA);
	Bloom::FileHeader headerB = Bloom::readHeader(tB);
	// Need to assert after every read operation
	assert(tA);
	assert(tB);

	const size_t IO_BUFFER_SIZE = 32 * 1024;
	unsigned char mask = 1;
	// The number of total bits in the vector
	size_t bitsA = headerA.endBitPos - headerA.startBitPos + 1;
	size_t bitsB = headerB.endBitPos - headerB.startBitPos + 1;
	// They need to be the same size to be comparable
	if(bitsA != bitsB ) {
		std::cerr << "Bit sizes of arrays not equal" << std::endl;
		exit(EXIT_FAILURE);
	}
	if (opt::verbose)
	std::cerr << "Bits: " << bitsA << std::endl;
	/* As in Choi et al. (2010),
	 a - cases where both bits are set (1/1)
	 b - cases where bits are set in the first but nor the second (1/0)
	 c - cases where bits are set in the second but not the first (0/1)
	 d - cases where bits are not set in either (0/0)
	 */
	unsigned long a = 0;
	unsigned long b = 0;
	unsigned long c = 0;
	unsigned long d = 0;
	// Iteratively compare bits
	for(size_t i = 0; i < bitsA;){
		char bufferA[IO_BUFFER_SIZE];
		char bufferB[IO_BUFFER_SIZE];
		// The number of bits in the buffer is its size * 8 except for the last iteration
		size_t bitsRead = std::min(IO_BUFFER_SIZE * 8, bitsA - i);
		size_t bytesRead = (bitsRead + 7)/8;
		// Read bytes from the the istream and immediately assert
		tA.read(bufferA, bytesRead);
		tB.read(bufferB, bytesRead);
		assert(tA);
		assert(tB);
		// For each byte in the buffer, compare bits
		for(size_t j = 0; j < IO_BUFFER_SIZE; j++){
			// Compare bit-wise
			for(int bit = 0; bit < 8; bit++){
				bool f = (bufferA[j] & (mask << bit)) != 0;
				bool s = (bufferB[j] & (mask << bit)) != 0;
				if( f == 1 && s == 1 ) {
					a++;
				} else if( f == 1 && s == 0) {
					b++;
				} else if( f == 0 && s == 1) {
					c++;
				} else d++;
			}
		}
		i += bitsRead;
	}
	assert(tA);
	assert(tB);
	// Result output:
	std::cout << "1/1: " << a << "\n1/0: " << b << "\n0/1: " << c << "\n0/0: " << d << std::endl;
	if(method == "jaccard"){
		float Dist = (float)a/(float)(a+b+c);
		std::cout << "Jaccard similarity: " << Dist << std::endl;
	}
	if(method == "czekanowski"){
		float Dist = (2*(float)a)/(float)((2*a)+b+c);
		std::cout << "Czekanowski similarity: " << Dist << std::endl;
	}
	if(method == "forbes"){
		float n = (float)(a + b + c + d);
		float Dist = (n*a - ((a+b)*(a+c))) / (n*std::min(a+b,a+c) - ((a+b) * (a+c)));
		std::cout << "Forbes similarity: " << Dist << std::endl;
	}
	// Check and clean up
	assert_good(tA, pathA);
	assert_good(tA, pathB);
	closeInputStream(inA, pathA);
	closeInputStream(inB, pathB);

  return 1;
}

int memberOf(int argc, char ** argv){
	// Initalise bloom and get globals
	BloomFilter bloom;
	parseGlobalOpts(argc, argv);
	// Arg parser to get `m' option in case set
	for (int c; (c = getopt_long(argc, argv,
								 shortopts, longopts, NULL)) != -1;) {
		istringstream arg(optarg != NULL ? optarg : "");
		switch (c) {
			case '?':
				cerr << PROGRAM ": unrecognized option: `-" << optopt
					<< "'" << endl;
				dieWithUsageError();
			case 'r':
				opt::inverse = true; break;
				break;
			case OPT_BED:
				opt::format = BED;
				break;
			case OPT_FASTA:
				opt::format = FASTA;
				break;
			case OPT_RAW:
				opt::format = RAW;
				break;
		}
		if (optarg != NULL && (!arg.eof() || arg.fail())) {
			cerr << PROGRAM ": invalid option: `-"
			<< (char)c << optarg << "'\n";
			exit(EXIT_FAILURE);
		}
	}
	string path = argv[optind];
	string fasta = argv[++optind];
	unsigned k = opt::k;
	if (opt::verbose)
		std::cerr << "Loading bloom filter from `"
		<< path << "'...\n";

	istream* in = openInputStream(path);
	assert_good(*in, path);
	*in >> bloom;

	assert(!fasta.empty());
	if (opt::verbose)
		std::cerr << "Reading `" << fasta << "'...\n";
	FastaReader _in(fasta.c_str(), FastaReader::FOLD_CASE);

	size_t seqCount=0;
	for (FastaRecord rec; _in >> rec; ++seqCount) {
		string& seq = rec.seq;
		if (seq.size() < k)
			continue;
		for (size_t i = 0; i < seq.size() - k + 1; ++i) {
			string kmer = seq.substr(i, k);
			size_t pos = kmer.find_last_not_of("ACGTacgt");
			if (pos != string::npos) {
				i += pos;
				continue;
			}
			if (bloom[Kmer(kmer)] || opt::inverse) {
				if (opt::format == FASTA) {
					cout << ">" << rec.id << ":seq:" << seqCount
						<< ":kmer:" << i << "\n";
				} else if (opt::format == BED) {
					cout << rec.id
						<< "\t" << i
						<< "\t" << i + k - 1
						<< "\t";
				}
				cout << kmer << "\n";
			}
		}
		if (opt::verbose && seqCount % 1000 == 0)
			cerr << "processed " << seqCount << " sequences" << endl;
	}
	assert(_in.eof());
	if (opt::verbose)
		cerr << "processed " << seqCount << " sequences" << endl;

	return 0;
}

/**
 * Calculate number of bases to trim from left end of sequence.
 */
int calcLeftTrim(const Sequence& seq, unsigned k, const BloomFilter& bloom,
	size_t minBranchLen)
{
	// Boost graph interface for Bloom filter
	DBGBloom<BloomFilter> g(bloom);

	// if this is the first k-mer we have found in
	// Bloom filter, starting from the left end
	// of the sequence
	bool firstKmerMatch = true;

	KmerIterator it(seq, k);
	for (; it != KmerIterator::end(); ++it) {

		const Kmer& kmer = *it;

		// assume k-mers not present in Bloom filter are
		// due to sequencing errors and should be trimmed
		if (!bloom[kmer])
			continue;

		// in degree, disregarding false branches
		unsigned inDegree = trueBranches(kmer, REVERSE, g,
				minBranchLen).size();
		// out degree, disregarding false branches
		unsigned outDegree = trueBranches(kmer, FORWARD, g,
				minBranchLen).size();

		if (firstKmerMatch) {
			bool leftTip = (inDegree == 0 && outDegree == 1);
			bool rightTip = (inDegree == 1 && outDegree == 0);
			if (!leftTip && !rightTip)
				break;
		} else if (inDegree != 1 || outDegree != 1) {
			// end of linear path
			break;
		}

		firstKmerMatch = false;

	} // for each k-mer (left to right)

	if (it.pos() == 0)
		return 0;

	return k + it.pos() - 1;
}

/**
 * Trim reads that corresponds to tips in the Bloom filter
 * de Bruijn graph.
 */
int trim(int argc, char** argv)
{
	// parse command line opts
	parseGlobalOpts(argc, argv);
	unsigned k = opt::k;

	// arg 1: Bloom filter
	// args 2-n: FASTA/FASTQ files
	if (argc - optind < 2) {
		cerr << PROGRAM ": missing arguments\n";
		dieWithUsageError();
	}

	// load Bloom filter de Bruijn graph
	string bloomPath(argv[optind++]);
	if (opt::verbose)
		cerr << "Loading bloom filter from `"
			<< bloomPath << "'...\n";

	BloomFilter bloom;
	istream *in = openInputStream(bloomPath);
	assert_good(*in, bloomPath);
	bloom.read(*in);
	assert_good(*in, bloomPath);

	if (opt::verbose)
		printBloomStats(cerr, bloom);

	// Calculate min length threshold for a "true branch"
	// (not due to Bloom filter false positives)
	const double falseBranchProbability = 0.0001;
	const size_t minBranchLen =
		(size_t)ceil(log(falseBranchProbability)/log(bloom.FPR()));

	if (opt::verbose >= 2)
		cerr << "min length threshold for true branches (k-mers): "
			<< minBranchLen << endl;

	// trim reads and print to STDOUT
	for (int i = optind; i < argc; ++i) {

		FastaReader in(argv[i], FastaReader::FOLD_CASE);
		for (FastqRecord rec; in >> rec;) {

			Sequence& seq = rec.seq;
			string& qual = rec.qual;

			// can't trim if read length < k; just echo
			// back to STDOUT
			if (seq.size() < k) {
				cout << rec;
				continue;
			}

			// start pos for trimmed read
			unsigned startPos = calcLeftTrim(seq, k, bloom, minBranchLen);
			// end pos for trimmed read
			unsigned endPos = seq.length() - 1 -
				calcLeftTrim(reverseComplement(seq), k, bloom, minBranchLen);

			// if whole read was trimmed away
			if (endPos < startPos)
				continue;

			// output trimmed read
			unsigned trimmedLen = endPos - startPos + 1;
			seq = seq.substr(startPos, trimmedLen);
			qual = qual.substr(startPos, trimmedLen);
			cout << rec;

		} // for each read
		assert(in.eof());

	} // for each input FASTA/FASTQ file

	// success
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
		exit(EXIT_SUCCESS);
	}
	if (command == "--version") {
		cout << VERSION_MESSAGE;
		exit(EXIT_SUCCESS);
	}
	else if (command == "build") {
		return build(argc, argv);
	}
	else if (command == "union") {
		return combine(argc, argv, BITWISE_OR);
	}
	else if (command == "intersect") {
		return combine(argc, argv, BITWISE_AND);
	}
	else if (command == "info") {
		return info(argc, argv);
	}
	else if (command == "compare") {
		return compare(argc, argv);
	}
	else if (command == "kmers" || command == "getKmers") {
		return memberOf(argc, argv);
	}
	else if (command == "trim") {
		return trim(argc, argv);
	}

	cerr << PROGRAM ": unrecognized command: `" << command
		<< "'" << endl;
	dieWithUsageError();
}
