#include "config.h"

#include "BloomDBG/bloom-dbg.h"
#include "BloomDBG/AssemblyCounters.h"
#include "BloomDBG/AssemblyParams.h"
#include "BloomDBG/Checkpoint.h"
#include "BloomDBG/HashAgnosticCascadingBloom.h"
#include "BloomDBG/MaskedKmer.h"
#include "BloomDBG/SpacedSeed.h"
#include "Common/StringUtil.h"
#include "Common/Options.h"
#include "DataLayer/Options.h"
#include "lib/bloomfilter/BloomFilter.hpp"

#include <getopt.h>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <iomanip>
#include <cstring>
#include <limits>
#include <string>

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
"    [options] <FASTQ> [FASTQ]... > assembly.fasta\n"
"\n"
"Perform a de Bruijn graph assembly of the given FASTQ files.\n"
"\n"
"Basic Options:\n"
"\n"
"  -b  --bloom-size=N           overall memory budget for the assembly in bytes.\n"
"                               Unit suffixes 'k' (kilobytes), 'M' (megabytes),\n"
"                               or 'G' (gigabytes) may be used. [required]\n"
"      --chastity               discard unchaste reads [default]\n"
"      --no-chastity            do not discard unchaste reads\n"
"  -g  --graph=FILE             write de Bruijn graph to FILE (GraphViz)\n"
"      --help                   display this help and exit\n"
"  -H  --num-hashes=N           number of Bloom filter hash functions [1]\n"
"  -i  --input-bloom=FILE       load Bloom filter from FILE\n"
"  -j, --threads=N              use N parallel threads [1]\n"
"      --trim-masked            trim masked bases from the ends of reads\n"
"      --no-trim-masked         do not trim masked bases from the ends\n"
"                               of reads [default]\n"
"  -k, --kmer=N                 the size of a k-mer [required]\n"
"      --kc=N                   use a cascading Bloom filter with N levels,\n"
"                               instead of a counting Bloom filter [2]\n"
"  -o, --out=FILE               write the contigs to FILE [STDOUT]\n"
"  -q, --trim-quality=N         trim bases from the ends of reads whose\n"
"                               quality is less than the threshold\n"
"  -Q, --mask-quality=N         mask all low quality bases as `N'\n"
"      --standard-quality       zero quality is `!' (33), typically\n"
"                               for FASTQ and SAM files [default]\n"
"      --illumina-quality       zero quality is `@' (64), typically\n"
"                               for qseq and export files\n"
"  -t, --trim-length            max branch length to trim, in k-mers [k]\n"
"  -v, --verbose                display verbose output\n"
"      --version                output version information and exit\n"
"\n"
"Spaced Seed Options:\n"
"\n"
"  -K, --single-kmer=N        use a spaced seed that consists of two k-mers\n"
"                             separated by a gap. K must be chosen such that\n"
"                             K <= k/2\n"
"      --qr-seed=N            use a spaced seed than consists of two mirrored\n"
"                             QR seeds separated by a gap.  The following must\n"
"                             hold: (a) N must be prime, (b) N >= 11,\n"
"                             (c) N <= k/2\n"
"  -s, --spaced-seed=STR      bitmask indicating k-mer positions to be\n"
"                             ignored during hashing. The pattern must be\n"
"                             symmetric\n"
"\n"
"Debugging Options:\n"
"\n"
"  -C, --cov-track=FILE       WIG track with 0/1 indicating k-mers with\n"
"                             coverage above the -c threshold. A reference\n"
"                             must also be specified with -R.\n"
"  -T, --trace-file=FILE      write debugging info about extension of\n"
"                             each read to FILE\n"
"  -R, --ref=FILE             specify a reference genome. FILE may be\n"
"                             FASTA, FASTQ, SAM, or BAM and may be gzipped."
"\n"
"Experimental Options:\n"
"\n"
"  Note!: These options may not be supported in future versions.\n"
"\n"
"      --checkpoint=N           create a checkpoint every N reads [disabled=0]\n"
"      --keep-checkpoint        do not delete checkpoint files after assembly\n"
"                               completes successfully [disabled]\n"
"      --checkpoint-prefix=STR  filename prefix for checkpoint files\n"
"                               ['bloom-dbg-checkpoint']\n"
"\n"
"Example:\n"
"\n"
"  Assemble a genome using a k-mer size of 50bp. Allocate a 1GB\n"
"  Bloom filter with 2 hash functions and require that a k-mer\n"
"  occurs 3 times or more to be included in the assembly. (The k-mer\n"
"  count threshold filters out k-mers containing sequencing errors.)\n"
"\n"
"  $ " PROGRAM " -k50 -b1G -H2 --kc=3 reads1.fq.gz reads2.fq.gz > assembly.fa\n"
"\n"
"Report bugs to <" PACKAGE_BUGREPORT ">.\n";

/** Assembly params (stores command-line options) */
BloomDBG::AssemblyParams params;

static const char shortopts[] = "b:C:g:H:i:j:k:K:o:q:Q:R:s:t:T:v";

enum {
	OPT_HELP = 1, OPT_VERSION, QR_SEED, MIN_KMER_COV,
	CHECKPOINT, KEEP_CHECKPOINT, CHECKPOINT_PREFIX
};

static const struct option longopts[] = {
	{ "bloom-size",        required_argument, NULL, 'b' },
	{ "min-coverage",      required_argument, NULL, 'c' },
	{ "cov-track",         required_argument, NULL, 'C' },
	{ "chastity",          no_argument, &opt::chastityFilter, 1 },
	{ "no-chastity",       no_argument, &opt::chastityFilter, 0 },
	{ "checkpoint",        required_argument, NULL, CHECKPOINT },
	{ "keep-checkpoint",   no_argument, NULL, KEEP_CHECKPOINT },
	{ "checkpoint-prefix", required_argument, NULL, CHECKPOINT_PREFIX },
	{ "graph",             required_argument, NULL, 'g' },
	{ "num-hashes",        required_argument, NULL, 'H' },
	{ "input-bloom",       required_argument, NULL, 'i' },
	{ "help",              no_argument, NULL, OPT_HELP },
	{ "threads",           required_argument, NULL, 'j' },
	{ "trim-masked",       no_argument, &opt::trimMasked, 1 },
	{ "no-trim-masked",    no_argument, &opt::trimMasked, 0 },
	{ "kmer",              required_argument, NULL, 'k' },
	{ "kc",                required_argument, NULL, MIN_KMER_COV },
	{ "single-kmer",       required_argument, NULL, 'K' },
	{ "out",               required_argument, NULL, 'o' },
	{ "trim-quality",      required_argument, NULL, 'q' },
	{ "mask-quality",      required_argument, NULL, 'Q' },
	{ "standard-quality",  no_argument, &opt::qualityOffset, 33 },
	{ "illumina-quality",  no_argument, &opt::qualityOffset, 64 },
	{ "qr-seed",           required_argument, NULL, QR_SEED },
	{ "ref",               required_argument, NULL, 'R' },
	{ "spaced-seed",       no_argument, NULL, 's' },
	{ "trim-length",       no_argument, NULL, 't' },
	{ "trace-file",        no_argument, NULL, 'T'},
	{ "verbose",           no_argument, NULL, 'v' },
	{ "version",           no_argument, NULL, OPT_VERSION },
	{ NULL, 0, NULL, 0 }
};

void printCascadingBloomStats(HashAgnosticCascadingBloom& bloom, ostream& os)
{
	for (unsigned i = 0; i < bloom.levels(); i++) {
		os << "Stats for Bloom filter level " << i+1 << ":\n"
			<< "\tBloom size (bits): "
			<< bloom.getBloomFilter(i).getFilterSize() << "\n"
			<< "\tBloom popcount (bits): "
			<< bloom.getBloomFilter(i).getPop() << "\n"
			<< "\tBloom filter FPR: " << setprecision(3)
			<< 100 * bloom.getBloomFilter(i).getFPR() << "%\n";
	}
}

/** Create optional auxiliary output files */
template <typename BloomFilterT>
void writeAuxiliaryFiles(int argc, char** argv, const BloomFilterT& bloom,
	const BloomDBG::AssemblyParams& params)
{
	/* generate wiggle coverage track */
	if (!params.covTrackPath.empty() && !params.refPath.empty())
		BloomDBG::writeCovTrack(bloom, params);

	/* generate de Bruijn graph in GraphViz format */
	if (!params.graphPath.empty()) {
		ofstream graphOut(params.graphPath.c_str());
		assert_good(graphOut, params.graphPath);
		BloomDBG::outputGraph(argc, argv, bloom, params, graphOut);
		assert_good(graphOut, params.graphPath);
		graphOut.close();
		assert_good(graphOut, params.graphPath);
	}
}

/** Initialize global variables for k-mer size and spaced seed pattern */
void initGlobals(const BloomDBG::AssemblyParams& params)
{
	/* set global variable for k-mer length */
	MaskedKmer::setLength(params.k);
	if (params.verbose)
		cerr << "Assembling with k-mer size " << params.k << endl;

	/* set global variable for spaced seed */
	if (params.K > 0)
		MaskedKmer::setMask(SpacedSeed::kmerPair(params.k, params.K));
	else if (params.qrSeedLen > 0)
		MaskedKmer::setMask(SpacedSeed::qrSeedPair(params.k, params.qrSeedLen));
	else
		MaskedKmer::setMask(params.spacedSeed);

	if (params.verbose && !MaskedKmer::mask().empty())
		cerr << "Using spaced seed " << MaskedKmer::mask() << endl;
}

/**
 * Resume assembly from previously saved checkpoint.
 */
void resumeAssemblyFromCheckpoint(int argc, char** argv,
	BloomDBG::AssemblyParams& params, ostream& out)
{
	assert(params.checkpointsEnabled() && checkpointExists(params));

	assert(params.initialized());
	initGlobals(params);

	/* empty Bloom filter de Bruijn graph */
	HashAgnosticCascadingBloom solidKmerSet;

	/* empty visited k-mers Bloom filter */
	BloomFilter visitedKmerSet;

	/* counters for progress messages */
	BloomDBG::AssemblyCounters counters;

	/* setup input/output streams for the assembly */

	/* input reads */
	FastaConcat in(argv + optind, argv + argc, FastaReader::FOLD_CASE);

	/* output stream for duplicate contigs FASTA output (for checkpoints) */
	ofstream checkpointOut;
	assert(!params.checkpointPathPrefix.empty());
	string prefix = params.checkpointPathPrefix;
	string checkpointPath = prefix + BloomDBG::CHECKPOINT_FASTA_EXT;
    checkpointOut.open(checkpointPath.c_str(), std::ofstream::app);
	assert_good(checkpointOut, checkpointPath);

	/* stream for trace file output ('-T' option) */
	ofstream traceOut;
	if (!params.tracePath.empty()) {
		traceOut.open(params.tracePath.c_str());
		assert_good(traceOut, params.tracePath);
		BloomDBG::SeqExtensionResult::printHeaders(traceOut);
		assert_good(traceOut, params.tracePath);
	}

	/* bundle input/output streams for assembly */
	BloomDBG::AssemblyStreams<FastaConcat>
		streams(in, out, checkpointOut, traceOut);

	/* restore state of Bloom filters, counters, and input/output streams */
	BloomDBG::resumeFromCheckpoint(solidKmerSet, visitedKmerSet,
		counters, params, streams);

	/* resume the assembly */
	BloomDBG::assemble(solidKmerSet, visitedKmerSet, counters, params,
		streams);
}

/**
 * Do the assembly after loading a pre-built Bloom filter
 * from file (`-i` option). (The input Bloom filter file
 * is constructed using `abyss-bloom build -t rolling-hash`.)
 */
void prebuiltBloomAssembly(int argc, char** argv,
	BloomDBG::AssemblyParams& params, ostream& out)
{
	/* load prebuilt Bloom filter from file */

	assert(!params.bloomPath.empty());
	if (params.verbose)
		cerr << "Loading prebuilt Bloom filter from `" << params.bloomPath
			<< "'" << endl;

	/* load the Bloom filter from file */

	HashAgnosticCascadingBloom bloom(params.bloomPath);

	if (params.verbose)
		cerr << "Bloom filter FPR: " << setprecision(3)
			<< bloom.FPR() * 100 << "%" << endl;

	printCascadingBloomStats(bloom, cerr);

	/* override command line options with values from Bloom file */

	params.k = bloom.getKmerSize();
	params.numHashes = bloom.getHashNum();
	params.bloomSize = bloom.size() + 7 / 8;
	if (params.trim == std::numeric_limits<unsigned>::max())
		params.trim = params.k;

	assert(params.numHashes <= MAX_HASHES);
	assert(params.initialized());

	/* init global vars for k-mer size and spaced seed pattern */

	initGlobals(params);

	if (params.verbose)
		cerr << params;

	/* do assembly */

	BloomDBG::assemble(argc - optind, argv + optind,
		bloom, params, out);

	/* write supplementary files (e.g. GraphViz) */

	writeAuxiliaryFiles(argc - optind, argv + optind, bloom, params);
}

/**
 * Load the reads into a cascading Bloom filter and do the assembly.
 */
void cascadingBloomAssembly(int argc, char** argv,
	const BloomDBG::AssemblyParams& params, ostream& out)
{
	/* init global vars for k-mer size and spaced seed pattern */

	assert(params.initialized());
	initGlobals(params);

	if (params.verbose)
		cerr << params;

	/*
	 * Build/load cascading Bloom filetr.
	 *
	 * Note 1: We use (params.minCov + 1) here because we use an additional
	 * Bloom filter in BloomDBG::assemble() to track the set of
	 * assembled k-mers.
	 *
	 * Note 2: BloomFilter class requires size to be a multiple of 64.
	 */

	const size_t bitsPerByte = 8;
	size_t bloomLevelSize = BloomDBG::roundUpToMultiple(
		params.bloomSize * bitsPerByte / (params.minCov + 1), (size_t)64);

	HashAgnosticCascadingBloom cascadingBloom(
		bloomLevelSize, params.numHashes, params.minCov, params.k);

	BloomDBG::loadBloomFilter(argc, argv, cascadingBloom, params.verbose);

	if (params.verbose)
		printCascadingBloomStats(cascadingBloom, cerr);

	/* second pass through FASTA files for assembling */

	BloomDBG::assemble(argc - optind, argv + optind,
		cascadingBloom, params, out);

	/* write supplementary files (e.g. GraphViz) */

	writeAuxiliaryFiles(argc - optind, argv + optind, cascadingBloom, params);
}

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
		  case 'C':
			arg >> params.covTrackPath; break;
		  case 'g':
			arg >> params.graphPath; break;
		  case 'H':
			arg >> params.numHashes; break;
		  case 'i':
			arg >> params.bloomPath; break;
		  case 'j':
			arg >> params.threads; break;
		  case 'k':
			arg >> params.k; break;
		  case 'K':
			params.resetSpacedSeedParams();
			arg >> params.K;
			break;
		  case 'o':
			arg >> params.outputPath; break;
		  case 'q':
			arg >> opt::qualityThreshold; break;
		  case 'R':
			arg >> params.refPath; break;
		  case 's':
			params.resetSpacedSeedParams();
			arg >> params.spacedSeed;
			break;
		  case 't':
			arg >> params.trim; break;
		  case 'T':
			arg >> params.tracePath; break;
		  case 'Q':
			arg >> opt::internalQThreshold; break;
		  case 'v':
			++params.verbose; break;
		  case OPT_HELP:
			cout << USAGE_MESSAGE;
			exit(EXIT_SUCCESS);
		  case MIN_KMER_COV:
			arg >> params.minCov; break;
		  case OPT_VERSION:
			cout << VERSION_MESSAGE;
			exit(EXIT_SUCCESS);
		  case QR_SEED:
			params.resetSpacedSeedParams();
			arg >> params.qrSeedLen;
			break;
		  case CHECKPOINT:
			arg >> params.readsPerCheckpoint;
			break;
		  case KEEP_CHECKPOINT:
			params.keepCheckpoint = true;
			break;
		  case CHECKPOINT_PREFIX:
			arg >> params.checkpointPathPrefix;
			break;
		}

		if (optarg != NULL && (!arg.eof() || arg.fail())) {
			cerr << PROGRAM ": invalid option: `-"
				<< (char)c << optarg << "'\n";
			exit(EXIT_FAILURE);
		}
	}

	if (params.bloomPath.empty() && params.bloomSize == 0) {
		cerr << PROGRAM ": missing mandatory option `-b'\n";
		die = true;
	}

	if (params.bloomPath.empty() && params.k == 0) {
		cerr << PROGRAM ": missing mandatory option `-k'\n";
		die = true;
	}

	if (params.k > 0 && params.K > 0 && params.K > params.k/2) {
		cerr << PROGRAM ": value of `-K' must be <= k/2\n";
		die = true;
	}

	if (params.numHashes > MAX_HASHES) {
		cerr << PROGRAM ": number of hash functions (`-H`) must "
			"be <= " << MAX_HASHES << " (set by `configure` option "
			"--enable-max-hashes=N)\n";
		die = true;
	}

	if (params.k > 0 && params.qrSeedLen > 0 &&
		(params.qrSeedLen < 11 || params.qrSeedLen > params.k/2)) {
		cerr << PROGRAM ": value of `--qr-seed' must be >= 11 and <= k/2\n";
		die = true;
	}

	if (!params.covTrackPath.empty() && params.refPath.empty()) {
		cerr << PROGRAM ": you must specify a reference with `-R' "
			"when using `-C'\n";
		die = true;
	}

	if (params.k > 0 && params.trim == std::numeric_limits<unsigned>::max()) {
		params.trim = params.k;
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

	/* print contigs to STDOUT unless -o option was set */
	ofstream outputFile;
	if (!params.outputPath.empty()) {
		outputFile.open(params.outputPath.c_str());
		assert_good(outputFile, params.outputPath);
	}
	ostream& out = params.outputPath.empty() ? cout : outputFile;

	/* load the Bloom filter and do the assembly */
	if (params.checkpointsEnabled() && checkpointExists(params))
		resumeAssemblyFromCheckpoint(argc, argv, params, out);
	else if (!params.bloomPath.empty())
		prebuiltBloomAssembly(argc, argv, params, out);
	else
		cascadingBloomAssembly(argc, argv, params, out);

	/* cleanup */
	if (!params.outputPath.empty())
		outputFile.close();

	return EXIT_SUCCESS;
}
