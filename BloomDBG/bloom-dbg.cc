#include "config.h"

#include "BloomDBG/bloom-dbg.h"
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
"    -G <genome_size> [options] <FASTQ> [FASTQ]... > assembly.fasta\n"
"\n"
"Perform a de Bruijn graph assembly of the given FASTQ files.\n"
"\n"
"Basic Options:\n"
"\n"
"  -b  --bloom-size=N         overall memory budget for the assembly in bytes.\n"
"                             Unit suffixes 'k' (kilobytes), 'M' (megabytes),\n"
"                             or 'G' (gigabytes) may be used. [required]\n"
"      --chastity             discard unchaste reads [default]\n"
"      --no-chastity          do not discard unchaste reads\n"
"  -g  --graph=FILE           write de Bruijn graph to FILE (GraphViz)\n"
"      --help                 display this help and exit\n"
"  -H  --num-hashes=N         number of Bloom filter hash functions [1]\n"
"  -j, --threads=N            use N parallel threads [1]\n"
"      --trim-masked          trim masked bases from the ends of reads\n"
"      --no-trim-masked       do not trim masked bases from the ends\n"
"                             of reads [default]\n"
"  -k, --kmer=N               the size of a k-mer [required]\n"
"      --kc=N                 use a cascading Bloom filter with N levels,\n"
"                             instead of a counting Bloom filter [2]\n"
"  -o, --out=FILE             write the contigs to FILE [STDOUT]\n"
"  -q, --trim-quality=N       trim bases from the ends of reads whose\n"
"                             quality is less than the threshold\n"
"  -Q, --mask-quality=N       mask all low quality bases as `N'\n"
"      --standard-quality     zero quality is `!' (33), typically\n"
"                             for FASTQ and SAM files [default]\n"
"      --illumina-quality     zero quality is `@' (64), typically\n"
"                             for qseq and export files\n"
"  -t, --trim-length          max branch length to trim, in k-mers [k]\n"
"  -v, --verbose              display verbose output\n"
"      --version              output version information and exit\n"
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

static const char shortopts[] = "b:C:g:H:j:k:K:o:q:Q:R:s:t:T:v";

enum { OPT_HELP = 1, OPT_VERSION, QR_SEED, MIN_KMER_COV };

static const struct option longopts[] = {
	{ "bloom-size",       required_argument, NULL, 'b' },
	{ "min-coverage",     required_argument, NULL, 'c' },
	{ "cov-track",        required_argument, NULL, 'C' },
	{ "chastity",         no_argument, &opt::chastityFilter, 1 },
	{ "no-chastity",      no_argument, &opt::chastityFilter, 0 },
	{ "graph",            required_argument, NULL, 'g' },
	{ "num-hashes",       required_argument, NULL, 'H' },
	{ "help",             no_argument, NULL, OPT_HELP },
	{ "threads",          required_argument, NULL, 'j' },
	{ "trim-masked",      no_argument, &opt::trimMasked, 1 },
	{ "no-trim-masked",   no_argument, &opt::trimMasked, 0 },
	{ "kmer",             required_argument, NULL, 'k' },
	{ "kc",               required_argument, NULL, MIN_KMER_COV },
	{ "single-kmer",      required_argument, NULL, 'K' },
	{ "out",              required_argument, NULL, 'o' },
	{ "trim-quality",     required_argument, NULL, 'q' },
	{ "mask-quality",     required_argument, NULL, 'Q' },
	{ "standard-quality", no_argument, &opt::qualityOffset, 33 },
	{ "illumina-quality", no_argument, &opt::qualityOffset, 64 },
	{ "qr-seed",          required_argument, NULL, QR_SEED },
	{ "ref",              required_argument, NULL, 'R' },
	{ "spaced-seed",      no_argument, NULL, 's' },
	{ "trim-length",      no_argument, NULL, 't' },
	{ "trace-file",       no_argument, NULL, 'T'},
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
			params.bloomSize = SIToBytes(arg); break;
		  case 'C':
			arg >> params.covTrackPath; break;
		  case 'g':
			arg >> params.graphPath; break;
		  case 'H':
			arg >> params.numHashes; break;
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
		}
		if (optarg != NULL && (!arg.eof() || arg.fail())) {
			cerr << PROGRAM ": invalid option: `-"
				<< (char)c << optarg << "'\n";
			exit(EXIT_FAILURE);
		}
	}

	if (params.bloomSize == 0) {
		cerr << PROGRAM ": missing mandatory option `-b'\n";
		die = true;
	}

	if (params.k == 0) {
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
		cerr << PROGRAM ": you must specify a reference with `-r' "
			"when using `-C'\n";
		die = true;
	}

	if (params.trim == std::numeric_limits<unsigned>::max()) {
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

	assert(params.initialized());

#if _OPENMP
	if (params.threads > 0)
		omp_set_num_threads(params.threads);
#endif

	/* set global variable for k-mer length */
	MaskedKmer::setLength(params.k);

	/* set global variable for spaced seed */
	if (params.K > 0)
		MaskedKmer::setMask(SpacedSeed::kmerPair(params.k, params.K));
	else if (params.qrSeedLen > 0)
		MaskedKmer::setMask(SpacedSeed::qrSeedPair(params.k, params.qrSeedLen));
	else
		MaskedKmer::setMask(params.spacedSeed);

	if (params.verbose && !MaskedKmer::mask().empty())
		cerr << "Using spaced seed " << MaskedKmer::mask() << endl;

	/* print contigs to STDOUT unless -o option was set */
	ofstream outputFile;
	if (!params.outputPath.empty()) {
		outputFile.open(params.outputPath.c_str());
		assert_good(outputFile, params.outputPath);
	}
	ostream& out = params.outputPath.empty() ? cout : outputFile;

	/* BloomFilter class requires size to be a multiple of 64 */
	const size_t bitsPerByte = 8;
	/*
	 * Note: it is (params.minCov + 1) here because we use an additional
	 * Bloom filter in BloomDBG::assemble() to track the set of
	 * assembled k-mers.
	 */
	size_t bloomLevelSize = BloomDBG::roundUpToMultiple(
		params.bloomSize * bitsPerByte / (params.minCov + 1), (size_t)64);

	/* use cascading Bloom filter to remove error k-mers */
	HashAgnosticCascadingBloom cascadingBloom(
		bloomLevelSize, params.numHashes, params.minCov, params.k);

	/* load reads into Bloom filter */
	for (int i = optind; i < argc; ++i) {
		/*
		 * Debugging feature: If there is a ':'
		 * separating the list of input read files into
		 * two parts, use the first set of files
		 * to load the Bloom filter and the second
		 * set of files for the assembly (read extension).
		 */
		if (strcmp(argv[i],":") == 0) {
			optind = i + 1;
			break;
		}
		BloomDBG::loadFile(cascadingBloom, argv[i], params.verbose);
	}
	if (params.verbose)
		cerr << "Bloom filter FPR: " << setprecision(3)
			<< cascadingBloom.FPR() * 100 << "%" << endl;

	if (!params.covTrackPath.empty()) {
		assert(!params.refPath.empty());
		BloomDBG::writeCovTrack(cascadingBloom, params);
	}

	/* second pass through FASTA files for assembling */
	BloomDBG::assemble(argc - optind, argv + optind,
		cascadingBloom, params, out);

	/* generate de Bruijn graph in GraphViz format (optional) */
	if (!params.graphPath.empty()) {
		ofstream graphOut(params.graphPath.c_str());
		assert_good(graphOut, params.graphPath);
		BloomDBG::outputGraph(argc - optind, argv + optind,
			cascadingBloom, params, graphOut);
		assert_good(graphOut, params.graphPath);
		graphOut.close();
		assert_good(graphOut, params.graphPath);
	}

	/* cleanup */
	if (!params.outputPath.empty())
		outputFile.close();

	return EXIT_SUCCESS;
}
