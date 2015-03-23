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
"  -P, --max-paths=N          merge at most N alternate paths; use 'nolimit'\n"
"                             for no limit [2]\n"
"  -q, --trim-quality=N       trim bases from the ends of reads whose\n"
"                             quality is less than the threshold\n"
"      --standard-quality     zero quality is `!' (33)\n"
"                             default for FASTQ and SAM files\n"
"      --illumina-quality     zero quality is `@' (64)\n"
"                             default for qseq and export files\n"
"  -r, --read-name=STR        only process reads with names that contain STR\n"
"  -s, --search-mem=N         mem limit for graph searches; multiply by the\n"
"                             number of threads (-j) to get the total mem used\n"
"                             for graph traversal [500M]\n"
"  -t, --trace-file=FILE      write graph search stats to FILE\n"
"  -v, --verbose              display verbose output\n"
"      --help                 display this help and exit\n"
"      --version              output version information and exit\n"
"\n"
"Report bugs to <" PACKAGE_BUGREPORT ">.\n";

const unsigned g_progressStep = 1000;
/*
 * ignore branches less than this length
 *(false positive branches)
 */
const unsigned g_trimLen = 3;

/*
 * Bloom filter use to keep track of portions
 * of genome that have already been assembled.
 * This Bloom filter is only used when the
 * -E (--extend) option is in effect.
 */
BloomFilter g_dupBloom;

namespace opt {

	/** The number of parallel threads. */
	static unsigned threads = 1;

	/** The size of the bloom filter in bytes. */
	size_t bloomSize = 500 * 1024 * 1024;

	/** The maximum count value of the BLoom filter. */
	unsigned max_count = 2;

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
	size_t mergedAndExtended;
	size_t mergedAndSkipped;
	size_t singleEndCorrected;
} g_count;

static const char shortopts[] = "b:B:d:D:eEf:F:i:Ij:k:lm:M:no:P:q:r:s:t:v";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
	{ "bloom-size",       required_argument, NULL, 'b' },
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
	{ "max-paths",        required_argument, NULL, 'P' },
	{ "trim-quality",     required_argument, NULL, 'q' },
	{ "standard-quality", no_argument, &opt::qualityOffset, 33 },
	{ "illumina-quality", no_argument, &opt::qualityOffset, 64 },
	{ "read-name",        required_argument, NULL, 'r' },
	{ "search-mem",       required_argument, NULL, 's' },
	{ "trace-file",       required_argument, NULL, 't' },
	{ "verbose",          no_argument, NULL, 'v' },
	{ "help",             no_argument, NULL, OPT_HELP },
	{ "version",          no_argument, NULL, OPT_VERSION },
	{ NULL, 0, NULL, 0 }
};

/**
 * Return true if the Bloom filter contains all of the
 * kmers in the given sequence.
 */
static bool bloomContainsSeq(const BloomFilter& bloom, const Sequence& seq)
{
	if (containsAmbiguityCodes(seq)) {
		Sequence seqCopy = seq;
		flattenAmbiguityCodes(seqCopy, false);
		for (KmerIterator it(seqCopy, opt::k); it != KmerIterator::end();
			++it) {
			if (!bloom[*it])
				return false;
		}
		return true;
	}
	for (KmerIterator it(seq, opt::k); it != KmerIterator::end(); ++it) {
		if (!bloom[*it])
			return false;
	}
	return true;
}

/**
 * Load the kmers of a given sequence into a Bloom filter.
 */
static inline void loadSeq(BloomFilter& bloom, unsigned k, const Sequence& seq)
{
	if (containsAmbiguityCodes(seq)) {
		Sequence seqCopy = seq;
		Sequence rc = reverseComplement(seqCopy);
		flattenAmbiguityCodes(seqCopy, false);
		flattenAmbiguityCodes(rc, false);
		Bloom::loadSeq(bloom, k, seqCopy);
		Bloom::loadSeq(bloom, k, rc);
	} else {
		Bloom::loadSeq(bloom, k, seq);
	}
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
static bool extendRead(Sequence& seq, unsigned k, const Graph& g)
{
	ExtendSeqResult result;
	bool extended = false;

	/*
	 * offset start pos to reduce chance of hitting
	 * a dead end on a false positive kmer
	 */
	const unsigned runLengthHint = 3;
	unsigned startPos = getStartKmerPos2(seq, k, FORWARD, g,
		runLengthHint);
	if (startPos != NO_MATCH) {
		assert(startPos <= seq.length() - k);
		result = extendSeq(seq, FORWARD, startPos, k, g,
				NO_LIMIT, g_trimLen, opt::mask);
		if (result == ES_EXTENDED_TO_DEAD_END ||
				result == ES_EXTENDED_TO_BRANCHING_POINT ||
				result == ES_EXTENDED_TO_CYCLE) {
			extended = true;
		}
	}

	startPos = getStartKmerPos2(seq, k, REVERSE, g, runLengthHint);
	if (startPos != NO_MATCH) {
		assert(startPos <= seq.length() - k);
		result = extendSeq(seq, REVERSE, startPos, k, g,
				NO_LIMIT, g_trimLen, opt::mask);
		if (result == ES_EXTENDED_TO_DEAD_END ||
				result == ES_EXTENDED_TO_BRANCHING_POINT ||
				result == ES_EXTENDED_TO_CYCLE) {
			extended = true;
		}
	}

	return extended;
}

enum ExtendResult { ER_NOT_EXTENDED, ER_REDUNDANT, ER_EXTENDED };

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
template <typename Graph>
static inline ExtendResult
extendReadIfNonRedundant(Sequence& seq, unsigned k, const Graph& g)
{
	bool redundant = false;
	if (opt::dupBloomSize > 0) {
		/*
		 * Check to see if the current pseudoread
		 * is contained in a region of the genome
		 * that has already been assembled.
		 */
#pragma omp critical(dupBloom)
		redundant = bloomContainsSeq(g_dupBloom, seq);
		if (redundant)
			return ER_REDUNDANT;
	}
	Sequence origSeq = seq;
	bool extended = extendRead(seq, k, g);
	if (opt::dupBloomSize > 0) {
		/*
		 * mark the extended read as an assembled
		 * region of the genome.
		 */
#pragma omp critical(dupBloom)
		{
			/* must check again to avoid race conditions */
			if (!bloomContainsSeq(g_dupBloom, origSeq))
				loadSeq(g_dupBloom, opt::k, seq);
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

/**
 * Print progress stats about reads merged/extended so far.
 */
static inline void printProgressMessage()
{
	cerr << "Merged " << g_count.uniquePath + g_count.multiplePaths << " of "
		<< g_count.readPairsProcessed << " read pairs";

	if (opt::extend) {
		cerr << ", corrected/extended " << g_count.singleEndCorrected << " of "
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


/**
 * For a successfully merged read pair, get the sequence
 * representing the connecting path between the two reads.
 */
template <typename Bloom>
static inline string getConnectingSeq(ConnectPairsResult& result,
	unsigned k, const Bloom& bloom)
{
	assert(result.pathResult == FOUND_PATH);
	(void)bloom;

	vector<FastaRecord>& paths = result.mergedSeqs;
	assert(paths.size() > 0);

	Sequence& seq = (paths.size() == 1) ?
		paths.front().seq : result.consensusSeq.seq;

	/*
	 * initialize sequence to the chars between the
	 * start and goal kmers of the path search.
	 */
	int startPos = result.startKmerPos;
	int endPos = seq.length() - result.goalKmerPos - k;
	assert(startPos >= 0 && startPos <=
		(int)(seq.length() - k + 1));

	return seq.substr(startPos, endPos - startPos + k);
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

	ConnectPairsResult result =
		connectPairs(opt::k, read1, read2, g, params);

	vector<FastaRecord>& paths = result.mergedSeqs;
	bool mergedSeqRedundant = false;
	bool read1Corrected = false;
	bool read1Redundant = false;
	bool read2Corrected = false;
	bool read2Redundant = false;

	/*
	 * extend reads inwards or outwards up to the
	 * next dead end or branching point in the de
	 * Brujin graph
	 */
	if (opt::extend) {
		ExtendResult extendResult;
		if (result.pathResult == FOUND_PATH
			&& result.pathMismatches <= params.maxPathMismatches
			&& result.readMismatches <= params.maxReadMismatches) {
			assert(paths.size() > 0);
			Sequence& seq = (paths.size() == 1) ?
				paths.front().seq : result.consensusSeq.seq;
			seq = getConnectingSeq(result, opt::k, bloom);
			extendResult = extendReadIfNonRedundant(
				seq, opt::k, g);
			if (extendResult == ER_REDUNDANT) {
#pragma omp atomic
				g_count.mergedAndSkipped++;
				mergedSeqRedundant = true;
			} else if (extendResult == ER_EXTENDED) {
#pragma omp atomic
				g_count.mergedAndExtended++;
			}
		} else {

			/*
			 * read pair could not be merged, so try
			 * to extend each read individually (in
			 * both directions).
			 */

//std::cerr << "correcting " << read1.id << " (read 1)" << std::endl;
			read1Corrected = correctAndExtendSeq(read1.seq,
				opt::k, g, read1.seq.length(), g_trimLen,
				opt::mask);

			if (read1Corrected) {
#pragma omp atomic
				g_count.singleEndCorrected++;
				extendResult = extendReadIfNonRedundant(read1.seq,
					opt::k, g);
				if (extendResult == ER_REDUNDANT)
					read1Redundant = true;
			}

//std::cerr << "correcting " << read2.id << " (read 2)" << std::endl;
			read2Corrected = correctAndExtendSeq(read2.seq,
				opt::k, g, read2.seq.length(), g_trimLen,
				opt::mask);

			if (read2Corrected) {
#pragma omp atomic
				g_count.singleEndCorrected++;
				extendResult = extendReadIfNonRedundant(read2.seq,
					opt::k, g);
				if (extendResult == ER_REDUNDANT)
					read2Redundant = true;
			}
		}
	}

	if (!opt::tracefilePath.empty())
#pragma omp critical(tracefile)
	{
		traceStream << result;
		assert_good(traceStream, opt::tracefilePath);
	}

	switch (result.pathResult) {

		case NO_PATH:
			assert(paths.empty());
			if (result.foundStartKmer && result.foundGoalKmer)
#pragma omp atomic
				++g_count.noPath;
			else {
#pragma omp atomic
				++g_count.noStartOrGoalKmer;
			}
			break;

		case FOUND_PATH:
			assert(!paths.empty());
			if (result.pathMismatches > params.maxPathMismatches ||
					result.readMismatches > params.maxReadMismatches) {
				if (result.pathMismatches > params.maxPathMismatches)
#pragma omp atomic
					++g_count.tooManyMismatches;
				else
					++g_count.tooManyReadMismatches;
				if (opt::extend) {
					if (read1Corrected || read2Corrected)
#pragma omp critical(mergedStream)
					{
						if (read1Corrected && !read1Redundant)
							mergedStream << (FastaRecord)read1;
						if (read2Corrected && !read2Redundant)
							mergedStream << (FastaRecord)read2;
					}
					if (!read1Corrected || !read2Corrected)
#pragma omp critical(readStream)
					{
						if (!read1Corrected)
							read1Stream << (FastaRecord)read1;
						if (!read2Corrected)
							read1Stream << (FastaRecord)read2;
					}
				} else
#pragma omp critical(readStream)
				{
					read1Stream << read1;
					read2Stream << read2;
				}
			}
			else if (paths.size() > 1) {
#pragma omp atomic
				++g_count.multiplePaths;
				if (!mergedSeqRedundant)
#pragma omp critical(mergedStream)
					mergedStream << result.consensusSeq;
			}
			else {
#pragma omp atomic
				++g_count.uniquePath;
				if (!mergedSeqRedundant)
#pragma omp critical(mergedStream)
					mergedStream << paths.front();
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

	if (result.pathResult != FOUND_PATH) {
		if (opt::extend) {
			if (read1Corrected || read2Corrected)
#pragma omp critical(mergedStream)
			{
				if (read1Corrected && !read1Redundant)
					mergedStream << (FastaRecord)read1;
				if (read2Corrected && !read2Redundant)
					mergedStream << (FastaRecord)read2;
			}
			if (!read1Corrected || !read2Corrected)
#pragma omp critical(readStream)
			{
				if (!read1Corrected)
					read1Stream << (FastaRecord)read1;
				if (!read2Corrected)
					read1Stream << (FastaRecord)read2;
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

	for (int c; (c = getopt_long(argc, argv,
					shortopts, longopts, NULL)) != -1;) {
		istringstream arg(optarg != NULL ? optarg : "");
		switch (c) {
		  case '?':
			die = true; break;
		  case 'b':
			opt::bloomSize = SIToBytes(arg); break;
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
		  case 'P':
			setMaxOption(opt::maxPaths, arg); break;
		  case 'q':
			arg >> opt::qualityThreshold; break;
		  case 'r':
			arg >> opt::readName; break;
		  case 's':
			opt::searchMem = SIToBytes(arg); break;
		  case 't':
			arg >> opt::tracefilePath; break;
		  case 'v':
			opt::verbose++; break;
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

#if _OPENMP
	if (opt::threads > 0)
		omp_set_num_threads(opt::threads);
#endif

	Kmer::setLength(opt::k);

#if USESEQAN
	seqanTests();
#endif

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

		// Specify bloom filter size in bits. Divide by two
		// because counting bloom filter requires twice as
		// much space.
		size_t bits = opt::bloomSize * 8 / opt::max_count;
		cascadingBloom = new CascadingBloomFilter(bits, opt::max_count);
#ifdef _OPENMP
		ConcurrentBloomFilter<CascadingBloomFilter> cbf(*cascadingBloom, 1000);
		for (int i = optind; i < argc; i++)
			Bloom::loadFile(cbf, opt::k, string(argv[i]), opt::verbose);
#else
		for (int i = optind; i < argc; i++)
			Bloom::loadFile(*cascadingBloom, opt::k, string(argv[i]), opt::verbose);
#endif
		bloom = &cascadingBloom->getBloomFilter(opt::max_count - 1);
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
	if (opt::extend)
		mergedOutputPath.append("_pseudoreads.fa");
	else
		mergedOutputPath.append("_merged.fa");
	ofstream mergedStream(mergedOutputPath.c_str());
	assert_good(mergedStream, mergedOutputPath);

	/*
	 * read pairs that were not successfully connected,
	 * but may have been extended inwards and/or outwards
	 * (if -E option was used)
	 */

	string read1OutputPath(opt::outputPrefix);
	if (opt::extend)
		read1OutputPath.append("_reads_1.fa");
	else
		read1OutputPath.append("_reads_1.fq");
	ofstream read1Stream(read1OutputPath.c_str());
	assert_good(read1Stream, read1OutputPath);

	string read2OutputPath(opt::outputPrefix);
	if (opt::extend)
		read2OutputPath.append("_reads_2.fa");
	else
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
	params.maxReadMismatches = opt::maxReadMismatches;
	params.kmerMatchesThreshold = 3;
	params.fixErrors = opt::fixErrors;
	params.maskBases = opt::mask;
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
					<< g_count.singleEndCorrected
					<< " (" << setprecision(3) <<  (float)100
					* g_count.singleEndCorrected / ((g_count.readPairsProcessed -
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
