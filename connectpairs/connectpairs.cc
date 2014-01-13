/**
 * Connect pairs using a Bloom filter de Bruijn graph
 * Copyright 2013 Shaun Jackman
 */

#include "config.h"

#include "DBGBloom.h"
#include "DBGBloomAlgorithms.h"

#include "Align/alignGlobal.h"
#include "Common/IOUtil.h"
#include "Common/Options.h"
#include "Common/StringUtil.h"
#include "DataLayer/FastaInterleave.h"
#include "DataLayer/FastaConcat.h"
#include "DataLayer/Options.h"
#include "Graph/DotIO.h"
#include "Graph/Options.h"
#include "Graph/GraphUtil.h"

#include <cassert>
#include <getopt.h>
#include <iostream>
#include <boost/tuple/tuple.hpp>

#if _OPENMP
# include <omp.h>
#endif

#undef USESEQAN

#if USESEQAN
#include <seqan/align.h>
#include <seqan/sequence.h>
#include <seqan/align_split.h>
#endif

using namespace std;
#if USESEQAN
using namespace seqan;
#endif
using boost::tie;

#define PROGRAM "abyss-connectpairs"

static const char VERSION_MESSAGE[] =
PROGRAM " (" PACKAGE_NAME ") " VERSION "\n"
"Written by Shaun Jackman, Hamid Mohamadi, Anthony Raymond and\n"
"Ben Vandervalk.\n"
"\n"
"Copyright 2013 Canada's Michael Smith Genome Science Centre\n";

static const char USAGE_MESSAGE[] =
"Usage: " PROGRAM " [OPTION]... [READS1 READS2]...\n"
"Connect the pairs READS1 and READS2 and close the gap using\n"
"a Bloom filter de Bruijn graph.\n"
"\n"
" Options:\n"
"\n"
"  -j, --threads=N            use N parallel threads [1]\n"
"  -k, --kmer=N               the size of a k-mer\n"
"  -b, --bloom-size=N         size of bloom filter [500M]\n"
"  -B, --max-branches=N       max branches in de Bruijn graph traversal [350]\n"
"  -f, --min-frag=N           min fragment size in base pairs [0]\n"
"  -F, --max-frag=N           max fragment size in base pairs [1000]\n"
"  -i, --input-bloom=FILE     load bloom filter from FILE\n"
"  -I, --interleaved          input reads files are interleaved\n"
"      --chastity             discard unchaste reads [default]\n"
"      --no-chastity          do not discard unchaste reads\n"
"      --trim-masked          trim masked bases from the ends of reads\n"
"      --no-trim-masked       do not trim masked bases from the ends\n"
"                             of reads [default]\n"
"  -M, --max-mismatches       max mismatches between merged paths [2]\n"
"  -o, --output-prefix=FILE   prefix of output FASTA files [required]\n"
"  -P, --max-paths=N          merge at most N alternate paths [2]\n"
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

	/** The size of the bloom filter in bytes. */
	size_t bloomSize = 500 * 1024 * 1024;

	/** Input read files are interleaved? */
	bool interleaved = false;

	/** Max active branches during de Bruijn graph traversal */
	unsigned maxBranches = 350;

	/** The size of a k-mer. */
	unsigned k;

	/** The minimum fragment size */
	unsigned minFrag = 0;

	/** The maximum fragment size */
	unsigned maxFrag = 1000;

	/** Bloom filter input file */
	static string inputBloomPath;

	/** Max paths between read 1 and read 2 */
	unsigned maxPaths = 2;

	/** Prefix for output files */
	static string outputPrefix;

	/** Max mismatches allowed when building consensus seqs */
	unsigned maxMismatches = 2;

}

/** Counters */
static struct {
	size_t noStartOrGoalKmer;
	size_t noPath;
	size_t uniquePath;
	size_t multiplePaths;
	size_t tooManyPaths;
	size_t tooManyBranches;
	size_t tooManyMismatches;
	size_t readPairsProcessed;
	size_t readPairsMerged;
} g_count;

static const char shortopts[] = "b:B:f:F:i:Ij:k:M:o:P:q:v";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
	{ "bloom-size",       required_argument, NULL, 'b' },
	{ "max-branches",     required_argument, NULL, 'B' },
	{ "min-frag",         required_argument, NULL, 'f' },
	{ "max-frag",         required_argument, NULL, 'F' },
	{ "input-bloom",      required_argument, NULL, 'i' },
	{ "interleaved",      no_argument, NULL, 'I' },
	{ "threads",          required_argument, NULL, 'j' },
	{ "kmer",             required_argument, NULL, 'k' },
	{ "chastity",         no_argument, &opt::chastityFilter, 1 },
	{ "no-chastity",      no_argument, &opt::chastityFilter, 0 },
	{ "trim-masked",      no_argument, &opt::trimMasked, 1 },
	{ "no-trim-masked",   no_argument, &opt::trimMasked, 0 },
	{ "output-prefix",    required_argument, NULL, 'o' },
	{ "max-mismatches",   required_argument, NULL, 'M' },
	{ "max-paths",        required_argument, NULL, 'P' },
	{ "trim-quality",     required_argument, NULL, 'q' },
	{ "standard-quality", no_argument, &opt::qualityOffset, 33 },
	{ "illumina-quality", no_argument, &opt::qualityOffset, 64 },
	{ "verbose",          no_argument, NULL, 'v' },
	{ "help",             no_argument, NULL, OPT_HELP },
	{ "version",          no_argument, NULL, OPT_VERSION },
	{ NULL, 0, NULL, 0 }
};

#if USESEQAN
const string r1 =
"AGAATCAACCAACCGTTCAATGATATAATCAAGAGCGATATTGTAATCTTTGTTTCT";
const string r2 =
"CGACGTCCACCAATTCGTCCCTGTGCACGAGCAGTTTCCAGTCCAGCTTTTGTTCGT";
const string ins =
"AGAATCAACCAACCGTTCAATGATATAATCAAGAGCGATATTGTAATCTTTGTTTCTGTCACCCGGCCCCCACGACTCAAGGATTAGACCATAAACACCATCCTCTTCACCTATCGAACACTCAGCTTTCAGTTCAATTCCATTATTATCAAAAACATGCATAATATTAATCTTTAATCAATTTTTCACGACAATACTACTTTTATTGATAAAATTGCAACAAGTTGCTGTTGTTTTACTTTCTTTTGTACACAAAGTGTCTTTAACTTTATTTATCCCCTGCAGGAAACCTCTTATACAAAGTTGACACACCAACATCATAGATAATCGCCACCTTCTGGCGAGGAGTTCCTGCTGCAATTAATCGTCCAGCTTGTGCCCATTGTTCTGGTGTAAGTTTGGGACGACGTCCACCAATTCGTCCCTGTGCACGAGCAGTTTCCAGTCCAGCTTTTGTTCGT";

static void seqanTests()
{
	typedef String<Dna> DS;
	typedef Align<DS> Alignment;

    //DS seq1 = "TTGT";
    //DS seq2 = "TTAGT";
	DS ref = ins;
	DS seq1 = r1;
	DS seq2 = r2;

    Alignment align1;
	resize(rows(align1), 2);
	assignSource(row(align1, 0), ref);
	assignSource(row(align1, 1), seq1);
    Alignment align2;
	resize(rows(align2), 2);
	assignSource(row(align2, 0), ref);
	assignSource(row(align2, 1), seq2);

	Score<int> scoring(2, -2, -50, -100);

	cout << splitAlignment(align1, align2, scoring) << endl;
	cout << align1 << endl;
	cout << align2 << endl;

	cout << localAlignment(align1, scoring) << endl;
	cout << align1 << endl;

	cout << localAlignment(align2, scoring) << endl;
	cout << align2 << endl;
}
#endif

/** Connect a read pair. */
static void connectPair(const DBGBloom& g,
	const FastqRecord& read1, const FastqRecord& read2,
	ofstream& mergedStream, ofstream& read1Stream,
	ofstream& read2Stream)
{
	SearchResult result
		= connectPairs(opt::k, read1, read2, g,
				opt::maxPaths, opt::minFrag,
				opt::maxFrag, opt::maxBranches);

	vector<FastaRecord>& paths = result.mergedSeqs;

#pragma omp atomic
	++g_count.readPairsProcessed;

	if (opt::verbose >= 2)
#pragma omp critical(cerr)
	{
		if(g_count.readPairsProcessed % g_progressStep == 0) {
			cerr << "Merged " << g_count.uniquePath + g_count.multiplePaths << " of "
				<< g_count.readPairsProcessed << " read pairs "
				<< "(no start/goal kmer: " << g_count.noStartOrGoalKmer << ", "
				<< "no path: " << g_count.noPath << ", "
				<< "too many paths: " << g_count.tooManyPaths << ", "
				<< "too many branches: " << g_count.tooManyBranches << ", "
				<< "too many mismatches: " << g_count.tooManyMismatches
				<< ")\n";
		}
	}

	switch (result.pathResult) {
	  case NO_PATH:
		assert(paths.empty());
		if (result.foundStartKmer && result.foundGoalKmer)
#pragma omp atomic
			++g_count.noPath;
		else
#pragma omp atomic
			++g_count.noStartOrGoalKmer;
		break;
	  case FOUND_PATH:
		assert(!paths.empty());
		if (paths.size() == 1) {
#pragma omp atomic
			++g_count.uniquePath;
#pragma omp critical(mergedStream)
			mergedStream << paths.front();
		} else {
			NWAlignment aln;
			unsigned matches, size;
			tie(matches, size) = align(paths, aln);
			assert(size >= matches);
			if (size - matches <= opt::maxMismatches) {
				FastaRecord read = paths.front();
				read.seq = aln.match_align;
#pragma omp atomic
				++g_count.multiplePaths;
#pragma omp critical(mergedStream)
				mergedStream << read;
			} else {
#pragma omp atomic
				++g_count.tooManyMismatches;
#pragma omp critical(read1Stream)
				read1Stream << read1;
#pragma omp critical(read2Stream)
				read2Stream << read2;
			}
		}
		break;
	  case TOO_MANY_PATHS:
#pragma omp atomic
		++g_count.tooManyPaths;
		break;
	  case TOO_MANY_BRANCHES:
#pragma omp atomic
		++g_count.tooManyBranches;
		break;
	}

	if (result.pathResult != FOUND_PATH) {
#pragma omp critical(read1Stream)
		read1Stream << read1;
#pragma omp critical(read2Stream)
		read2Stream << read2;
	}
}

/** Connect read pairs. */
template <typename FastaStream>
static void connectPairs(const DBGBloom& g, FastaStream& in,
	ofstream& mergedStream, ofstream& read1Stream, ofstream& read2Stream)
{
#pragma omp parallel
	for (FastqRecord a, b;;) {
		bool good;
#pragma omp critical(in)
		good = in >> a >> b;
		if (good)
			connectPair(g, a, b, mergedStream, read1Stream, read2Stream);
		else
			break;
	}
}

/**
 * Connect pairs using a Bloom filter de Bruijn graph
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
		  case 'B':
			arg >> opt::maxBranches; break;
		  case 'f':
			arg >> opt::minFrag; break;
		  case 'F':
			arg >> opt::maxFrag; break;
		  case 'i':
			arg >> opt::inputBloomPath; break;
		  case 'I':
			opt::interleaved = true; break;
		  case 'j':
			arg >> opt::threads; break;
		  case 'k':
			arg >> opt::k; break;
		  case 'M':
			arg >> opt::maxMismatches; break;
		  case 'o':
			arg >> opt::outputPrefix; break;
		  case 'P':
			arg >> opt::maxPaths; break;
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

	if (opt::outputPrefix.empty()) {
		cerr << PROGRAM ": missing mandatory option `-o'\n";
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

#if _OPENMP
	if (opt::threads > 0)
		omp_set_num_threads(opt::threads);
#endif

	Kmer::setLength(opt::k);

#if USESEQAN
	seqanTests();
#endif

	assert(opt::bloomSize > 0);

	BloomFilterBase* bloom = NULL;

	if (!opt::inputBloomPath.empty()) {

		if (opt::verbose)
			std::cerr << "Loading bloom filter from `"
				<< opt::inputBloomPath << "'...\n";

		const char* inputPath = opt::inputBloomPath.c_str();
		ifstream inputBloom(inputPath, ios_base::in | ios_base::binary);
		assert_good(inputBloom, inputPath);
		BloomFilter* loadedBloom = new BloomFilter();
		inputBloom >> *loadedBloom;
		assert_good(inputBloom, inputPath);
		inputBloom.close();
		bloom = loadedBloom;

	} else {

		// Specify bloom filter size in bits. Divide by two
		// because counting bloom filter requires twice as
		// much space.
		size_t bits = opt::bloomSize * 8 / 2;
		bloom = new CountingBloomFilter(bits);
		for (int i = optind; i < argc; i++)
			bloom->loadFile(opt::k, string(argv[i]), opt::verbose);

	}

	if (opt::verbose)
		cerr << "Bloom filter FPR: " << setprecision(3)
			<< 100 * bloom->FPR() << "%\n";

	DBGBloom g(*bloom);

	string mergedOutputPath(opt::outputPrefix);
	mergedOutputPath.append("_merged.fa");
	ofstream mergedStream(mergedOutputPath.c_str());
	assert_good(mergedStream, mergedOutputPath);

	string read1OutputPath(opt::outputPrefix);
	read1OutputPath.append("_reads_1.fq");
	ofstream read1Stream(read1OutputPath.c_str());
	assert_good(read1Stream, read1OutputPath);

	string read2OutputPath(opt::outputPrefix);
	read2OutputPath.append("_reads_2.fq");
	ofstream read2Stream(read2OutputPath.c_str());
	assert_good(read2Stream, read2OutputPath);

	if (opt::verbose > 0)
		cerr << "Connecting read pairs\n";

	if (opt::interleaved) {
		FastaConcat in(argv + optind, argv + argc,
				FastaReader::FOLD_CASE);
		connectPairs(g, in, mergedStream, read1Stream, read2Stream);
		assert(in.eof());
	} else {
		FastaInterleave in(argv + optind, argv + argc,
				FastaReader::FOLD_CASE);
		connectPairs(g, in, mergedStream, read1Stream, read2Stream);
		assert(in.eof());
	}

	if (opt::verbose > 0) {
		cerr <<
			"Processed " << g_count.readPairsProcessed << " read pairs\n"
			"Merged (Unique path + Multiple paths): "
				<< g_count.uniquePath + g_count.multiplePaths
				<< " (" << setprecision(3) <<  (float)100
				    * (g_count.uniquePath + g_count.multiplePaths) /
				   g_count.readPairsProcessed
				<< "%)\n"
			"No start/goal kmer: " << g_count.noStartOrGoalKmer
				<< " (" << setprecision(3) << (float)100
					* g_count.noStartOrGoalKmer / g_count.readPairsProcessed
				<< "%)\n"
			"No path: " << g_count.noPath
				<< " (" << setprecision(3) << (float)100
					* g_count.noPath / g_count.readPairsProcessed
				<< "%)\n"
			"Unique path: " << g_count.uniquePath
				<< " (" << setprecision(3) << (float)100
					* g_count.uniquePath / g_count.readPairsProcessed
				<< "%)\n"
			"Multiple paths: " << g_count.multiplePaths
				<< " (" << setprecision(3) << (float)100
					* g_count.multiplePaths / g_count.readPairsProcessed
				<< "%)\n"
			"Too many paths: " << g_count.tooManyPaths
				<< " (" << setprecision(3) << (float)100
					* g_count.tooManyPaths / g_count.readPairsProcessed
				<< "%)\n"
			"Too many branches: " << g_count.tooManyBranches
				<< " (" << setprecision(3) << (float)100
					* g_count.tooManyBranches / g_count.readPairsProcessed
				<< "%)\n"
			"Too many mismatches: " << g_count.tooManyMismatches
				<< " (" << setprecision(3) << (float)100
					* g_count.tooManyMismatches / g_count.readPairsProcessed
				<< "%)\n"
			"Bloom filter FPR: " << setprecision(3) << 100 * bloom->FPR()
				<< "%\n";
	}

	assert_good(mergedStream, mergedOutputPath.c_str());
	mergedStream.close();
	assert_good(read1Stream, read1OutputPath.c_str());
	read1Stream.close();
	assert_good(read2Stream, read2OutputPath.c_str());
	read2Stream.close();

	delete bloom;

	return 0;
}
