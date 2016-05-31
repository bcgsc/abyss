/**
 * Connect pairs using a Bloom filter de Bruijn graph
 * Copyright 2013 Shaun Jackman
 */

#include "config.h"

#include "konnector.h"
#include "Bloom/CascadingBloomFilter.h"
#include "DBGBloom.h"
#include "DBGBloomAlgorithms.h"

#include "Align/alignGlobal.h"
#include "Common/IOUtil.h"
#include "Common/Options.h"
#include "Common/StringUtil.h"
#include "DataLayer/FastaConcat.h"
#include "DataLayer/FastaInterleave.h"
#include "DataLayer/Options.h"
#include "Graph/DotIO.h"
#include "Graph/Options.h"
#include "Graph/GraphUtil.h"

#include <cassert>
#include <getopt.h>
#include <iostream>
#include <cstring>
#include <algorithm>

#if _OPENMP
# include <omp.h>
# include "Bloom/ConcurrentBloomFilter.h"
#endif

using namespace std;

#define PROGRAM "konnector"

static const char VERSION_MESSAGE[] =
PROGRAM " (" PACKAGE_NAME ") " VERSION "\n"
"Written by Shaun Jackman, Hamid Mohamadi, Anthony Raymond, \n"
"Ben Vandervalk and Justin Chu.\n"
"\n"
"Copyright 2014 Canada's Michael Smith Genome Science Centre\n";

static const char USAGE_MESSAGE[] =
"Usage: " PROGRAM " -k <kmer_size> -o <output_prefix> [options]... <reads1> [reads2]...\n"
"Connect the pairs READS1 and READS2 and close the gap using\n"
"a Bloom filter de Bruijn graph.\n"
"\n"
" Options:\n"
"\n"
"  -j, --threads=N            use N parallel threads [1]\n"
"  -k, --kmer=N               the size of a k-mer\n"
"  -b, --bloom-size=N         size of bloom filter [500M]\n"
"  -c, --min-coverage=N       kmer coverage threshold for error correction [2].\n"
"                             This option specifies the number of levels in the\n"
"                             cascading Bloom filter; it has no effect if the Bloom\n"
"                             filter is loaded from an external file.\n"
"  -B, --max-branches=N       max branches in de Bruijn graph traversal;\n"
"                             use 'nolimit' for no limit [350]\n"
"  -d, --dot-file=FILE        write graph traversals to a DOT file\n"
"  -D, --dup-bloom-size=N     use an additional Bloom filter to avoid\n"
"                             assembling the same region of the genome\n"
"                             multiple times. This option is highly\n"
"                             recommended when the -E (--extend) option\n"
"                             and has no effect otherwise. As a rule of\n"
"                             thumb, the Bloom filter size should be\n"
"                             about twice the target genome size [disabled]\n"
"  -e, --fix-errors           find and fix single-base errors when reads\n"
"                             have no kmers in bloom filter [disabled]\n"
"  -E, --extend               in addition to finding a connecting path,\n"
"                             extend the reads outwards to the next\n"
"                             dead end or branching point in the de Brujin\n"
"                             graph. If the reads were not successfully\n"
"                             connected, extend them inwards as well.\n"
"  --fastq                    output merged reads in FASTQ format\n"
"                             (default is FASTA)\n"
"  -f, --min-frag=N           min fragment size in base pairs [0]\n"
"  -F, --max-frag=N           max fragment size in base pairs [1000]\n"
"  -i, --input-bloom=FILE     load bloom filter from FILE\n"
"  -I, --interleaved          input reads files are interleaved\n"
"      --mask                 mask new and changed bases as lower case\n"
"      --no-mask              do not mask bases [default]\n"
"      --chastity             discard unchaste reads [default]\n"
"      --no-chastity          do not discard unchaste reads\n"
"      --trim-masked          trim masked bases from the ends of reads\n"
"      --no-trim-masked       do not trim masked bases from the ends\n"
"                             of reads [default]\n"
"  -m, --read-mismatches=N    max mismatches between paths and reads; use\n"
"                             'nolimit' for no limit [nolimit]\n"
"  -M, --max-mismatches=N     max mismatches between all alternate paths;\n"
"                             use 'nolimit' for no limit [2]\n"
"  -n  --no-limits            disable all limits; equivalent to\n"
"                             '-B nolimit -m nolimit -M nolimit -P nolimit'\n"
"  -o, --output-prefix=FILE   prefix of output FASTA files [required]\n"
"  --preserve-reads           don't correct any bases within the reads [disabled]\n"
"  -p, --alt-paths-mode       output a separate pseudoread for each alternate\n"
"                             path connecting a read pair (default is to create\n"
"                             a consensus sequence of all connecting paths)\n"
"  -P, --max-paths=N          merge at most N alternate paths; use 'nolimit'\n"
"                             for no limit [2]\n"
"  -q, --trim-quality=N       trim bases from the ends of reads whose\n"
"                             quality is less than the threshold\n"
"      --standard-quality     zero quality is `!' (33), typically\n"
"                             for FASTQ and SAM files [default]\n"
"      --illumina-quality     zero quality is `@' (64), typically\n"
"                             for qseq and export files\n"
"  -Q, --corrected-qual       quality score for bases corrected or inserted\n"
"                             by konnector; only relevant when --fastq is\n"
"                             in effect [40]\n"
"  -r, --read-name=STR        only process reads with names that contain STR\n"
"  -s, --search-mem=N         mem limit for graph searches; multiply by the\n"
"                             number of threads (-j) to get the total mem used\n"
"                             for graph traversal [500M]\n"
"  -t, --trace-file=FILE      write graph search stats to FILE\n"
"  -v, --verbose              display verbose output\n"
"  -x, --read-identity=N      min percent seq identity between consensus seq\n"
"                             and reads [0]\n"
"  -X, --path-identity=N      min percent seq identity across alternate\n"
"                             connecting paths [0]\n"
"      --help                 display this help and exit\n"
"      --version              output version information and exit\n"
"\n"
"Report bugs to <" PACKAGE_BUGREPORT ">.\n";

const unsigned g_progressStep = 1000;
/**
 * ignore branches less than this length
 *(false positive branches)
 */
const unsigned g_trimLen = 3;
/*
 * Bloom filter to keep track of portions
 * of genome that have already been assembled.
 * This Bloom filter is only used when both
 * the --extend and --dup-bloom-size options
 * are in effect.
 */
BloomFilter g_dupBloom;

namespace opt {

	/** The number of parallel threads. */
	static unsigned threads = 1;

	/** The size of the bloom filter in bytes. */
	size_t bloomSize = 500 * 1024 * 1024;

	/** The maximum count value of the Bloom filter. */
	unsigned minCoverage = 2;

	/** Input read files are interleaved? */
	bool interleaved = false;

	/** Max active branches during de Bruijn graph traversal */
	unsigned maxBranches = 350;

	/** multi-graph DOT file containing graph traversals */
	static string dotPath;

	/**
	 * Dup Bloom filter size.
	 * The dup filter is used to avoid assembling duplicate
	 * sequences when the -E (--extend) option is in effect.
	 */
	size_t dupBloomSize = 0;

	/**
	 * Find and fix single base errors when a read has no
	 * kmers in the bloom filter.
	 */
	bool fixErrors = false;

	/**
	 * Extend reads outwards until the next dead or branching
	 * point in the de Bruijn graph.  If a read pair is not
	 * successfully connected, extend them inwards as well.
	 */
	bool extend = false;

	/**
	 * Output pseudo-reads in FASTQ format.
	 */
	bool fastq = false;

	/** The size of a k-mer. */
	unsigned k;

	/** The minimum fragment size */
	unsigned minFrag = 0;

	/** The maximum fragment size */
	unsigned maxFrag = 1000;

	/** Bloom filter input file */
	static string inputBloomPath;

	/**
	 * Do not correct bases in input read sequences.
	 */
	static bool preserveReads = false;

	/**
	 * Output separate sequence for each alternate path
	 * between read pairs
	 */
	static bool altPathsMode = false;

	/** Max paths between read 1 and read 2 */
	unsigned maxPaths = 2;

	/**
	 * Quality score for bases that are corrected
	 * or inserted by konnector.
	 */
	uint8_t correctedQual = 40;

	/** Prefix for output files */
	static string outputPrefix;

	/** Max mismatches allowed when building consensus seqs */
	unsigned maxMismatches = 2;

	/** Only process reads that contain this substring. */
	static string readName;

	/** Max mem used per thread during graph traversal */
	static size_t searchMem = 500 * 1024 * 1024;

	/** Output file for graph search stats */
	static string tracefilePath;

	/** Mask bases not in reads */
	static int mask = 0;

	/** Max mismatches between consensus and original reads */
	static unsigned maxReadMismatches = NO_LIMIT;

	/**
	 * Min percent seq identity between consensus seq
	 * and input reads
	 */
	static float minReadIdentity = 0.0f;

	/**
	 * Min percent seq identity between all alternate
	 * paths
	 */
	static float minPathIdentity = 0.0f;
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
	size_t tooManyReadMismatches;
	size_t containsCycle;
	size_t exceededMemLimit;
	size_t traversalMemExceeded;
	size_t readPairsProcessed;
	size_t readPairsMerged;
	size_t skipped;
	/* counts below are used only when -E is enabled */
	size_t mergedAndSkipped;
	size_t singleEndExtended;
} g_count;

static const char shortopts[] = "b:B:c:d:D:eEf:F:i:Ij:k:lm:M:no:p:P:q:Q:r:s:t:vx:X:";

enum { OPT_FASTQ = 1, OPT_HELP, OPT_PRESERVE_READS, OPT_VERSION };

static const struct option longopts[] = {
	{ "bloom-size",       required_argument, NULL, 'b' },
	{ "min-coverage",     required_argument, NULL, 'c' },
	{ "max-branches",     required_argument, NULL, 'B' },
	{ "dot-file",         required_argument, NULL, 'd' },
	{ "dup-bloom-size",   required_argument, NULL, 'D' },
	{ "fix-errors",       no_argument, NULL, 'e' },
	{ "extend",           no_argument, NULL, 'E' },
	{ "min-frag",         required_argument, NULL, 'f' },
	{ "max-frag",         required_argument, NULL, 'F' },
	{ "input-bloom",      required_argument, NULL, 'i' },
	{ "interleaved",      no_argument, NULL, 'I' },
	{ "threads",          required_argument, NULL, 'j' },
	{ "kmer",             required_argument, NULL, 'k' },
	{ "chastity",         no_argument, &opt::chastityFilter, 1 },
	{ "no-chastity",      no_argument, &opt::chastityFilter, 0 },
	{ "mask",             no_argument, &opt::mask, 1 },
	{ "no-mask",          no_argument, &opt::mask, 0 },
	{ "no-limits",        no_argument, NULL, 'n' },
	{ "trim-masked",      no_argument, &opt::trimMasked, 1 },
	{ "no-trim-masked",   no_argument, &opt::trimMasked, 0 },
	{ "output-prefix",    required_argument, NULL, 'o' },
	{ "read-mismatches",  required_argument, NULL, 'm' },
	{ "max-mismatches",   required_argument, NULL, 'M' },
	{ "alt-paths-mode",   no_argument, NULL, 'p' },
	{ "max-paths",        required_argument, NULL, 'P' },
	{ "trim-quality",     required_argument, NULL, 'q' },
	{ "corrected-qual",   required_argument, NULL, 'Q' },
	{ "standard-quality", no_argument, &opt::qualityOffset, 33 },
	{ "illumina-quality", no_argument, &opt::qualityOffset, 64 },
	{ "read-name",        required_argument, NULL, 'r' },
	{ "search-mem",       required_argument, NULL, 's' },
	{ "trace-file",       required_argument, NULL, 't' },
	{ "verbose",          no_argument, NULL, 'v' },
	{ "read-identity",    required_argument, NULL, 'x' },
	{ "path-identity",    required_argument, NULL, 'X' },
	{ "fastq",            no_argument, NULL, OPT_FASTQ },
	{ "help",             no_argument, NULL, OPT_HELP },
	{ "preserve-reads",   no_argument, NULL, OPT_PRESERVE_READS },
	{ "version",          no_argument, NULL, OPT_VERSION },
	{ NULL, 0, NULL, 0 }
};

/**
 * Return true if the Bloom filter contains all of the
 * "good" kmers in the given sequence.
 */
static inline bool isSeqRedundant(const BloomFilter& assembledKmers,
	const BloomFilter& goodKmers, Sequence seq)
{
	flattenAmbiguityCodes(seq, false);
	for (KmerIterator it(seq, opt::k); it != KmerIterator::end(); ++it) {
		if (goodKmers[*it] && !assembledKmers[*it])
			return false;
	}
	return true;
}

/**
 * Load the kmers of a given sequence into a Bloom filter.
 */
static inline void addKmers(BloomFilter& bloom,
	const BloomFilter& goodKmers, unsigned k,
	const Sequence& seq)
{
	if (containsAmbiguityCodes(seq)) {
		Sequence flattened = seq;
		Sequence rcFlattened = reverseComplement(seq);
		flattenAmbiguityCodes(flattened, false);
		flattenAmbiguityCodes(rcFlattened, false);
		for (KmerIterator it(flattened, k);
			it != KmerIterator::end();++it) {
			if (goodKmers[*it])
				bloom.insert(*it);
		}
		for (KmerIterator it(rcFlattened, k);
			it != KmerIterator::end(); ++it) {
			if (goodKmers[*it])
				bloom.insert(*it);
		}
		return;
	} else {
		for (KmerIterator it(seq, k);
			it != KmerIterator::end(); ++it) {
			if (goodKmers[*it])
				bloom.insert(*it);
		}
	}
}

enum ExtendResult { ER_NOT_EXTENDED, ER_REDUNDANT, ER_EXTENDED };

/**
 * Calculate quality string for a pseudo-read.  A base will
 * have a score of CORRECTED_BASE_QUAL if it was corrected
 * by konnector or added by konnector (in the gap between
 * paired-end reads).  For bases that are unchanged from the
 * input reads, the original quality score is used.  In the
 * case that the two input read(s) overlap and both provide
 * a correct base call, the maximum of the two quality scores
 * is used.
 */
static inline std::string calcQual(const FastqRecord& read1,
	const FastqRecord& read2, Sequence& merged)
{
	unsigned char correctedQual = opt::qualityOffset + opt::correctedQual;
	std::string qual(merged.length(), correctedQual);

	/*
	 * In the case that the input files are FASTA,
	 * the quality strings for read1 / read2 will be
	 * empty, so just return a uniform quality string.
	 */
	if (read1.qual.empty() || read2.qual.empty())
		return qual;

	Sequence r1 = read1.seq, r2 = reverseComplement(read2.seq);
	std::string r1qual = read1.qual, r2qual = read2.qual;
	std::reverse(r2qual.begin(), r2qual.end());
	assert(r1.length() <= merged.length());
	assert(r2.length() <= merged.length());

	/* region covered only by read 1 */
	unsigned r2offset = merged.length() - r2.length();
	for (unsigned r1pos = 0; r1pos < r1.length() && r1pos < r2offset;
		++r1pos) {
		if (r1.at(r1pos) == merged.at(r1pos)) {
			qual.at(r1pos) = r1qual.at(r1pos);
		} else {
			//r1Corrected.at(i) = true;
			qual.at(r1pos) = correctedQual;
		}
	}

	/* region where read 1 and read 2 overlap */
	for (unsigned r1pos = r2offset; r1pos < r1.length(); ++r1pos) {
		unsigned r2pos = r1pos - r2offset;
		if (r1.at(r1pos) != merged.at(r1pos) ||
			r2.at(r2pos) != merged.at(r1pos)) {
			qual.at(r1pos) = correctedQual;
		} else {
			assert(r1.at(r1pos) == r2.at(r2pos));
			qual.at(r1pos) = max(r1qual.at(r1pos), r2qual.at(r2pos));
		}
	}

	/* region covered only by read 2 */
	for (unsigned r1pos = max(r2offset, (unsigned)r1.length());
		r1pos < merged.length(); ++r1pos) {
		unsigned r2pos = r1pos - r2offset;
		if (r2.at(r2pos) == merged.at(r1pos)) {
			qual.at(r1pos) = r2qual.at(r2pos);
		} else {
			qual.at(r1pos) = correctedQual;
		}
	}

	return qual;
}

static inline string calcQual(const FastqRecord& orig,
	const Sequence& extended, unsigned extendedLeft,
	unsigned extendedRight)
{
	assert(extended.length() == orig.seq.length() +
		extendedLeft + extendedRight);

	unsigned char correctedQual = opt::qualityOffset + opt::correctedQual;
	string qual(extended.length(), correctedQual);

	/*
	 * In the case that the input files are FASTA,
	 * the quality strings for read1 / read2 will be
	 * empty, so just return a uniform quality string.
	 */
	if (orig.qual.empty())
		return qual;

	unsigned offset = extendedLeft;
	for (unsigned i = 0; i < orig.seq.length(); ++i) {
		assert(offset + i < extended.length());
		assert(i < orig.seq.length());
		assert(i < orig.qual.length());
		if (orig.seq.at(i) == extended.at(offset + i))
			qual.at(offset + i) = orig.qual.at(i);
	}

	return qual;
}

/**
 * Extend a read/pseudoread both left and right until
 * we hit the next dead end or branching point in the
 * de Bruijn graph.
 *
 * @param seq sequence to be extended
 * @param k kmer size
 * @param g de Bruijn graph
 * return true if the read was extended in either
 * (or both) directions, false otherwise
 */
template <typename Graph>
static bool extendRead(FastqRecord& rec, unsigned k, const Graph& g)
{
	unsigned extendedLeft = 0, extendedRight = 0;
	Sequence extendedSeq = rec.seq;

	/*
	 * offset start pos to reduce chance of hitting
	 * a dead end on a false positive kmer
	 */
	const unsigned runLengthHint = 3;
	unsigned startPos = getStartKmerPos(extendedSeq, k, FORWARD, g,
		runLengthHint);
	if (startPos != NO_MATCH) {
		assert(startPos <= extendedSeq.length() - k);
		unsigned lengthBefore = extendedSeq.length();
		extendSeq(extendedSeq, FORWARD, startPos, k, g,
			NO_LIMIT, g_trimLen, opt::mask,
			!opt::altPathsMode, opt::preserveReads);
		extendedRight = extendedSeq.length() - lengthBefore;
	}

	startPos = getStartKmerPos(extendedSeq, k, REVERSE, g, runLengthHint);
	if (startPos != NO_MATCH) {
		assert(startPos <= extendedSeq.length() - k);
		unsigned lengthBefore = extendedSeq.length();
		extendSeq(extendedSeq, REVERSE, startPos, k, g,
			NO_LIMIT, g_trimLen, opt::mask,
			!opt::altPathsMode, opt::preserveReads);
		extendedLeft = extendedSeq.length() - lengthBefore;
	}

	if (extendedLeft > 0 || extendedRight > 0) {
		rec.qual = calcQual(rec, extendedSeq,
			extendedLeft, extendedRight);
		rec.seq = extendedSeq;
		return true;
	}

	return false;
}

/**
 * Attempt to extend a merged read (a.k.a. pseudoread)
 * outward to the next branching point or dead end in
 * the de Bruijn graph.
 *
 * @param seq pseudoread to be extended
 * @param k kmer size
 * @param g de Bruijn graph in which to perform extension
 * @return ExtendResult (ER_NOT_EXTENDED, ER_EXTENDED,
 * ER_REDUNDANT)
 */
template <typename Graph, typename BloomT1, typename BloomT2>
static inline ExtendResult
extendReadIfNonRedundant(FastqRecord& seq, BloomT1& assembledKmers,
	const BloomT2& goodKmers, unsigned k, const Graph& g)
{
	bool extended = false;
	bool redundant = false;

	if (opt::dupBloomSize > 0) {
		/*
		 * Check to see if the current pseudoread
		 * is contained in a region of the genome
		 * that has already been assembled.
		 */
#pragma omp critical(dupBloom)
		redundant = isSeqRedundant(assembledKmers, goodKmers, seq);
		if (redundant)
			return ER_REDUNDANT;
	}
	Sequence origSeq = seq.seq;
	extended = extendRead(seq, k, g);
	if (opt::dupBloomSize > 0) {
		/*
		 * mark the extended read as an assembled
		 * region of the genome.
		 */
#pragma omp critical(dupBloom)
		{
			/* must check again to avoid race conditions */
			if (!isSeqRedundant(assembledKmers, goodKmers, origSeq))
				addKmers(assembledKmers, goodKmers, k, seq.seq);
			else
				redundant = true;
		}
		if (redundant)
			return ER_REDUNDANT;
	}
	assert(!redundant);
	if (extended)
		return ER_EXTENDED;
	else
		return ER_NOT_EXTENDED;
}

static inline FastqRecord connectingSeq(const FastqRecord& mergedSeq,
	unsigned startKmerPos, unsigned goalKmerPos)
{
	FastqRecord rec;

	unsigned start = startKmerPos;
	unsigned end = mergedSeq.seq.length() - 1 - goalKmerPos;
	assert(start <= end);

	rec.id = mergedSeq.id;
	rec.seq = mergedSeq.seq.substr(start, end - start + 1);
	rec.qual = mergedSeq.qual.substr(start, end - start + 1);

	return rec;
}

/**
 * Print progress stats about reads merged/extended so far.
 */
static inline void printProgressMessage()
{
	cerr << "Merged " << g_count.uniquePath + g_count.multiplePaths << " of "
		<< g_count.readPairsProcessed << " read pairs";

	if (opt::extend) {
		cerr << ", corrected/extended " << g_count.singleEndExtended << " of "
			<< (g_count.readPairsProcessed - g_count.uniquePath -
				g_count.multiplePaths) * 2
		<< " unmerged reads";
	}

	cerr << " (no start/goal kmer: " << g_count.noStartOrGoalKmer << ", "
		<< "no path: " << g_count.noPath << ", "
		<< "too many paths: " << g_count.tooManyPaths << ", "
		<< "too many branches: " << g_count.tooManyBranches << ", "
		<< "too many path/path mismatches: " << g_count.tooManyMismatches << ", "
		<< "too many path/read mismatches: " << g_count.tooManyReadMismatches << ", "
		<< "contains cycle: " << g_count.containsCycle << ", "
		<< "skipped: " << g_count.skipped
		<< ")\n";
}

static inline void updateCounters(const ConnectPairsParams& params,
	const ConnectPairsResult& result)
{
	switch (result.pathResult) {
		case NO_PATH:
			assert(result.mergedSeqs.empty());
			if (result.foundStartKmer && result.foundGoalKmer)
#pragma omp atomic
				++g_count.noPath;
			else
#pragma omp atomic
				++g_count.noStartOrGoalKmer;
			break;

		case FOUND_PATH:
			assert(!result.mergedSeqs.empty());
			if (result.pathMismatches > params.maxPathMismatches ||
				result.pathIdentity < params.minPathIdentity) {
#pragma omp atomic
					++g_count.tooManyMismatches;
			} else if (result.readMismatches > params.maxReadMismatches ||
				result.readIdentity < params.minReadIdentity) {
#pragma omp atomic
					++g_count.tooManyReadMismatches;
			} else {
				if (result.mergedSeqs.size() == 1)
#pragma omp atomic
					++g_count.uniquePath;
				else
#pragma omp atomic
					++g_count.multiplePaths;
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

		case PATH_CONTAINS_CYCLE:
#pragma omp atomic
			++g_count.containsCycle;
			break;

		case EXCEEDED_MEM_LIMIT:
#pragma omp atomic
			++g_count.exceededMemLimit;
			break;
	}
}

static inline void outputRead(const FastqRecord& read, ostream& out,
	bool fastq = true)
{
	if (fastq)
		out << read;
	else
		out << (FastaRecord)read;
}

static inline bool exceedsMismatchThresholds(const ConnectPairsParams& params,
	const ConnectPairsResult& result)
{
	return (result.pathMismatches > params.maxPathMismatches ||
		result.pathIdentity < params.minPathIdentity ||
		result.readMismatches > params.maxReadMismatches ||
		result.readIdentity < params.minReadIdentity);
}

/**
 * Correct and extend an unmerged single-end read.
 * @return true if the read was modified, false otherwise
 */
template <typename Graph, typename BloomT1, typename BloomT2>
static inline bool correctAndExtend(FastqRecord& read,
	BloomT1& assembledKmers, const BloomT2& goodKmers,
	unsigned k, const Graph& g, bool preserveRead=false)
{
	bool corrected = false;
	if (!preserveRead)
		corrected = trimRead(read, k, g);
	if (preserveRead || corrected) {
		ExtendResult extendResult =
			extendReadIfNonRedundant(read, assembledKmers,
				goodKmers, k, g);
		if (extendResult == ER_EXTENDED)
			return true;
	}
	return corrected;
}

/** Connect a read pair. */
template <typename Graph, typename Bloom>
static void connectPair(const Graph& g,
	const Bloom& bloom,
	FastqRecord& read1,
	FastqRecord& read2,
	const ConnectPairsParams& params,
	ofstream& mergedStream,
	ofstream& read1Stream,
	ofstream& read2Stream,
	ofstream& traceStream)
{
	/*
	 * Implements the -r option, which is used to only
	 * process a subset of the input read pairs.
	 */
	if (!opt::readName.empty() &&
		read1.id.find(opt::readName) == string::npos) {
#pragma omp atomic
		++g_count.skipped;
		return;
	}

	/* Search for connecting paths between read pair */

	ConnectPairsResult result =
		connectPairs(opt::k, read1, read2, g, params);

	/* Calculate quality strings for merged reads */

	vector<FastqRecord> paths;
	FastqRecord consensus;
	if (result.pathResult == FOUND_PATH) {
		for (unsigned i = 0; i < result.mergedSeqs.size(); ++i) {
			FastqRecord fastq;
			fastq.id = result.mergedSeqs.at(i).id;
			fastq.seq = result.mergedSeqs.at(i).seq;
			fastq.qual = calcQual(read1, read2, result.mergedSeqs.at(i).seq);
			paths.push_back(fastq);
		}
		consensus.id = result.consensusSeq.id;
		consensus.seq = result.consensusSeq.seq;
		consensus.qual = calcQual(read1, read2, result.consensusSeq.seq);
	}

	bool outputRead1 = false;
	bool outputRead2 = false;
	std::vector<bool> pathRedundant;

	/*
	 * extend reads inwards or outwards up to the
	 * next dead end or branching point in the de
	 * Brujin graph
	 */
	if (opt::extend) {
		ExtendResult extendResult;
		if (result.pathResult == FOUND_PATH &&
			!exceedsMismatchThresholds(params, result)) {
			/* we found at least one connecting path */
			assert(paths.size() > 0);
			if (opt::altPathsMode) {
				/* extend each alternate path independently */
				for (unsigned i = 0; i < paths.size(); ++i) {
					if (!opt::preserveReads)
						paths.at(i) = connectingSeq(paths.at(i),
							result.startKmerPos, result.goalKmerPos);
					extendResult = extendReadIfNonRedundant(
						paths.at(i), g_dupBloom, bloom, opt::k, g);
					pathRedundant.push_back(extendResult == ER_REDUNDANT);
				}
			} else  {
				/* extend consensus sequence for all paths */
				if (!opt::preserveReads)
					consensus = connectingSeq(consensus,
						result.startKmerPos, result.goalKmerPos);
				extendResult = extendReadIfNonRedundant(
					consensus, g_dupBloom, bloom, opt::k, g);
				pathRedundant.push_back(extendResult == ER_REDUNDANT);
			}
			if (std::find(pathRedundant.begin(), pathRedundant.end(),
				false) == pathRedundant.end()) {
#pragma omp atomic
				g_count.mergedAndSkipped++;
			}
		} else {

			/*
			 * read pair could not be merged, so try
			 * to correct and extend each read individually (in
			 * both directions).
			 */

			if (correctAndExtend(read1, g_dupBloom, bloom,
				opt::k, g, opt::preserveReads)) {
					/* avoid duplicate read IDs */
					if (!endsWith(read1.id, "/1")) {
						read1.id.append("/1");
						read1.comment.clear();
					}
					outputRead1 = true;
#pragma omp atomic
					g_count.singleEndExtended++;
			}

			if (correctAndExtend(read2, g_dupBloom, bloom,
				opt::k, g, opt::preserveReads)) {
					/* avoid duplicate read IDs */
					if (!endsWith(read2.id, "/2")) {
						read2.id.append("/2");
						read2.comment.clear();
					}
					outputRead2 = true;
#pragma omp atomic
					g_count.singleEndExtended++;
			}

		}
	}

	if (!opt::tracefilePath.empty())
#pragma omp critical(tracefile)
	{
		traceStream << result;
		assert_good(traceStream, opt::tracefilePath);
	}

	/* update stats regarding merge successes / failures */

	updateCounters(params, result);

	/* ouput merged / unmerged reads */

	if (result.pathResult == FOUND_PATH &&
		!exceedsMismatchThresholds(params, result)) {
		assert(!paths.empty());
		if (opt::altPathsMode) {
#pragma omp critical(mergedStream)
			for (unsigned i = 0; i < paths.size(); ++i) {
				if (opt::dupBloomSize == 0 || !pathRedundant.at(i))
					outputRead(paths.at(i), mergedStream, opt::fastq);
			}
		} else if (opt::dupBloomSize == 0 || !pathRedundant.front()) {
#pragma omp critical(mergedStream)
			outputRead(consensus, mergedStream, opt::fastq);
		}
	} else {
		if (opt::extend) {
			if (outputRead1 || outputRead2)
#pragma omp critical(mergedStream)
			{
				if (outputRead1)
					outputRead(read1, mergedStream, opt::fastq);
				if (outputRead2)
					outputRead(read2, mergedStream, opt::fastq);
			}
			if (!outputRead1 || !outputRead2)
#pragma omp critical(readStream)
			{
				if (!outputRead1)
					read1Stream << read1;
				if (!outputRead2)
					read2Stream << read2;
			}
		} else
#pragma omp critical(readStream)
		{
			read1Stream << read1;
			read2Stream << read2;
		}
	}
}

/** Connect read pairs. */
template <typename Graph, typename FastaStream, typename Bloom>
static void connectPairs(const Graph& g,
	const Bloom& bloom,
	FastaStream& in,
	const ConnectPairsParams& params,
	ofstream& mergedStream,
	ofstream& read1Stream,
	ofstream& read2Stream,
	ofstream& traceStream)
{
#pragma omp parallel
	for (FastqRecord a, b;;) {
		bool good;
#pragma omp critical(in)
		good = in >> a >> b;
		if (good) {
			connectPair(g, bloom, a, b, params, mergedStream, read1Stream,
				read2Stream, traceStream);
#pragma omp atomic
			g_count.readPairsProcessed++;
			if (opt::verbose >= 2)
#pragma omp critical(cerr)
			{
				if(g_count.readPairsProcessed % g_progressStep == 0)
					printProgressMessage();
			}
		} else {
			break;
		}
	}
}

/**
 * Set the value for a commandline option, using "nolimit"
 * to represent NO_LIMIT.
 */
static inline void setMaxOption(unsigned& arg, istream& in)
{
	string str;
	getline(in, str);
	if (!in.fail() && str.compare("nolimit")==0) {
		arg = NO_LIMIT;
	} else {
		istringstream ss(str);
		ss >> arg;
		// copy state bits (fail, bad, eof) to
		// original stream
		in.clear(ss.rdstate());
	}
}

/**
 * Connect pairs using a Bloom filter de Bruijn graph
 */
int main(int argc, char** argv)
{
	bool die = false;
	bool minCovOptUsed = false;

	for (int c; (c = getopt_long(argc, argv,
					shortopts, longopts, NULL)) != -1;) {
		istringstream arg(optarg != NULL ? optarg : "");
		switch (c) {
		  case '?':
			die = true; break;
		  case 'b':
			opt::bloomSize = SIToBytes(arg); break;
		  case 'c':
			arg >> opt::minCoverage;
			minCovOptUsed = true;
			break;
		  case 'B':
			setMaxOption(opt::maxBranches, arg); break;
		  case 'd':
			arg >> opt::dotPath; break;
		  case 'D':
			opt::dupBloomSize = SIToBytes(arg); break;
		  case 'e':
			opt::fixErrors = true; break;
		  case 'E':
			opt::extend = true; break;
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
		  case 'm':
			setMaxOption(opt::maxReadMismatches, arg); break;
		  case 'n':
			opt::maxBranches = NO_LIMIT;
			opt::maxReadMismatches = NO_LIMIT;
			opt::maxMismatches = NO_LIMIT;
			opt::maxPaths = NO_LIMIT;
			break;
		  case 'M':
			setMaxOption(opt::maxMismatches, arg); break;
		  case 'o':
			arg >> opt::outputPrefix; break;
		  case 'p':
			opt::altPathsMode = true; break;
		  case 'P':
			setMaxOption(opt::maxPaths, arg); break;
		  case 'q':
			arg >> opt::qualityThreshold; break;
		  case 'Q':
			arg >> opt::correctedQual; break;
		  case 'r':
			arg >> opt::readName; break;
		  case 's':
			opt::searchMem = SIToBytes(arg); break;
		  case 't':
			arg >> opt::tracefilePath; break;
		  case 'x':
			arg >> opt::minReadIdentity; break;
		  case 'X':
			arg >> opt::minPathIdentity; break;
		  case 'v':
			opt::verbose++; break;
		  case OPT_FASTQ:
			opt::fastq = true; break;
		  case OPT_HELP:
			cout << USAGE_MESSAGE;
			exit(EXIT_SUCCESS);
		  case OPT_PRESERVE_READS:
			opt::preserveReads = true; break;
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

	if (opt::outputPrefix.empty()) {
		cerr << PROGRAM ": missing mandatory option `-o'\n";
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

	if (!opt::inputBloomPath.empty() && minCovOptUsed) {
		cerr << PROGRAM ": warning: -c option has no effect when "
			" using a pre-built Bloom filter (-i option)\n";
	}


#if _OPENMP
	if (opt::threads > 0)
		omp_set_num_threads(opt::threads);
#endif

	Kmer::setLength(opt::k);

#if USESEQAN
	seqanTests();
#endif

	/*
	 * We need to set a default quality score offset
	 * in order to generate quality scores
	 * for bases that are corrected/inserted by
	 * konnector (--fastq option).
	 */
	if (opt::qualityOffset == 0)
		opt::qualityOffset = 33;

	assert(opt::bloomSize > 0);

	if (opt::dupBloomSize > 0)
		g_dupBloom.resize(opt::dupBloomSize * 8);

	BloomFilter* bloom;
	CascadingBloomFilter* cascadingBloom = NULL;

	if (!opt::inputBloomPath.empty()) {

		if (opt::verbose)
			std::cerr << "Loading bloom filter from `"
				<< opt::inputBloomPath << "'...\n";

		bloom = new BloomFilter();

		const char* inputPath = opt::inputBloomPath.c_str();
		ifstream inputBloom(inputPath, ios_base::in | ios_base::binary);
		assert_good(inputBloom, inputPath);
		inputBloom >> *bloom;
		assert_good(inputBloom, inputPath);
		inputBloom.close();

	} else {

		if (opt::verbose)
			std::cerr << "Using a minimum kmer coverage threshold of "
				<< opt::minCoverage << "\n";

		// Specify bloom filter size in bits and divide by number
		// of levels in cascading Bloom filter.

		size_t bits = opt::bloomSize * 8 / opt::minCoverage;
		cascadingBloom = new CascadingBloomFilter(bits, opt::minCoverage);
#ifdef _OPENMP
		ConcurrentBloomFilter<CascadingBloomFilter> cbf(*cascadingBloom, 1000);
		for (int i = optind; i < argc; i++)
			Bloom::loadFile(cbf, opt::k, string(argv[i]), opt::verbose);
#else
		for (int i = optind; i < argc; i++)
			Bloom::loadFile(*cascadingBloom, opt::k, string(argv[i]), opt::verbose);
#endif
		bloom = &cascadingBloom->getBloomFilter(opt::minCoverage - 1);
	}

	if (opt::verbose)
		cerr << "Bloom filter FPR: " << setprecision(3)
			<< 100 * bloom->FPR() << "%\n";

	ofstream dotStream;
	if (!opt::dotPath.empty()) {
		if (opt::verbose)
			cerr << "Writing graph traversals to "
				"dot file `" << opt::dotPath << "'\n";
		dotStream.open(opt::dotPath.c_str());
		assert_good(dotStream, opt::dotPath);
	}

	ofstream traceStream;
	if (!opt::tracefilePath.empty()) {
		if (opt::verbose)
			cerr << "Writing graph search stats to `"
				<< opt::tracefilePath << "'\n";
		traceStream.open(opt::tracefilePath.c_str());
		assert(traceStream.is_open());
		ConnectPairsResult::printHeaders(traceStream);
		assert_good(traceStream, opt::tracefilePath);
	}

	DBGBloom<BloomFilter> g(*bloom);

	/*
	 * read pairs that were successfully connected
	 * (and possibly extended outwards)
	 */

	string mergedOutputPath(opt::outputPrefix);
	mergedOutputPath.append("_pseudoreads");
	if (opt::fastq)
		mergedOutputPath.append(".fq");
	else
		mergedOutputPath.append(".fa");
	ofstream mergedStream(mergedOutputPath.c_str());
	assert_good(mergedStream, mergedOutputPath);

	/*
	 * read pairs that were not successfully connected,
	 * but may have been extended inwards and/or outwards
	 * (if -E option was used)
	 */

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

	ConnectPairsParams params;

	params.minMergedSeqLen = opt::minFrag;
	params.maxMergedSeqLen = opt::maxFrag;
	params.maxPaths = opt::maxPaths;
	params.maxBranches = opt::maxBranches;
	params.maxPathMismatches = opt::maxMismatches;
	params.minPathIdentity = opt::minPathIdentity;
	params.maxReadMismatches = opt::maxReadMismatches;
	params.minReadIdentity = opt::minReadIdentity;
	params.kmerMatchesThreshold = 3;
	params.fixErrors = opt::fixErrors;
	params.maskBases = opt::mask;
	params.preserveReads = opt::preserveReads;
	params.memLimit = opt::searchMem;
	params.dotPath = opt::dotPath;
	params.dotStream = opt::dotPath.empty() ? NULL : &dotStream;

	if (opt::interleaved) {
		FastaConcat in(argv + optind, argv + argc,
				FastaReader::FOLD_CASE);
		connectPairs(g, *bloom, in, params, mergedStream, read1Stream,
				read2Stream, traceStream);
		assert(in.eof());
	} else {
		FastaInterleave in(argv + optind, argv + argc,
				FastaReader::FOLD_CASE);
		connectPairs(g, *bloom, in, params, mergedStream, read1Stream,
				read2Stream, traceStream);
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
			"Too many path/path mismatches: " << g_count.tooManyMismatches
				<< " (" << setprecision(3) << (float)100
					* g_count.tooManyMismatches / g_count.readPairsProcessed
				<< "%)\n"
			"Too many path/read mismatches: " << g_count.tooManyReadMismatches
				<< " (" << setprecision(3) << (float)100
					* g_count.tooManyReadMismatches / g_count.readPairsProcessed
				<< "%)\n"
			"Contains cycle: " << g_count.containsCycle
				<< " (" << setprecision(3) << (float)100
					* g_count.containsCycle / g_count.readPairsProcessed
				<< "%)\n"
			"Skipped: " << g_count.skipped
				<< " (" << setprecision(3) << (float)100
					* g_count.skipped / g_count.readPairsProcessed
				<< "%)\n";
			if (opt::extend) {
				cerr << "Unmerged reads corrected/extended: "
					<< g_count.singleEndExtended
					<< " (" << setprecision(3) <<  (float)100
					* g_count.singleEndExtended / ((g_count.readPairsProcessed -
					g_count.uniquePath - g_count.multiplePaths) * 2)
					<< "%)\n";
			}
			std::cerr << "Bloom filter FPR: " << setprecision(3)
				<< 100 * bloom->FPR() << "%\n";
	}

	if (!opt::inputBloomPath.empty())
		delete bloom;
	else
		delete cascadingBloom;

	assert_good(mergedStream, mergedOutputPath.c_str());
	mergedStream.close();
	assert_good(read1Stream, read1OutputPath.c_str());
	read1Stream.close();
	assert_good(read2Stream, read2OutputPath.c_str());
	read2Stream.close();

	if (!opt::dotPath.empty()) {
		assert_good(dotStream, opt::dotPath);
		dotStream.close();
	}

	if (!opt::tracefilePath.empty()) {
		assert_good(traceStream, opt::tracefilePath);
		traceStream.close();
	}

	return 0;
}
