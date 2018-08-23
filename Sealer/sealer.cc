/**
 * Close intra-scaffold gaps
 * Copyright 2014 Canada's Michael Smith Genome Science Centre
 */

#include "config.h"

#include "Konnector/konnector.h"
#include "Konnector/DBGBloom.h"
#include "Konnector/DBGBloomAlgorithms.h"
#include "Bloom/CascadingBloomFilter.h"

#include "Align/alignGlobal.h"
#include "Common/IOUtil.h"
#include "Common/Options.h"
#include "Common/StringUtil.h"
#include "DataLayer/FastaConcat.h"
#include "DataLayer/Options.h"
#include "Graph/DotIO.h"
#include "Graph/Options.h"
#include "Graph/GraphUtil.h"

#include <cassert>
#include <getopt.h>
#include <iostream>
#include <cstring>

#include <string>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <map>


#if _OPENMP
# include <omp.h>
# include "Bloom/ConcurrentBloomFilter.h"
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

#define PROGRAM "abyss-sealer"

static const char VERSION_MESSAGE[] =
PROGRAM " (" PACKAGE_NAME ") " VERSION "\n"
"Written by Shaun Jackman, Hamid Mohamadi, Anthony Raymond,\n"
"Ben Vandervalk and Daniel Paulino\n"
"\n"
"Copyright 2014 Canada's Michael Smith Genome Science Centre\n";

static const char USAGE_MESSAGE[] =
"Usage: " PROGRAM "-b <Bloom filter size> -k <kmer size> -k <kmer size>... -o <output_prefix> -S <path to scaffold file> [options]... <reads1> [reads2]...\n"
"i.e. abyss-sealer -b20G -k90 -k80 -k70 -k60 -k50 -k40 -k30 -o test -S scaffold.fa read1.fa read2.fa\n\n"
"Close gaps by using left and right flanking sequences of gaps as 'reads' for Konnector\n"
"and performing multiple runs with each of the supplied K values.\n"
"\n"
" Options:\n"
"\n"
"      --print-flanks           outputs flank files\n"
"  -S, --input-scaffold=FILE    load scaffold from FILE\n"
"  -L, --flank-length=N         length of flanks to be used as pseudoreads [100]\n"
"  -D, --flank-distance=N       distance of flank from gap [0]\n"
"  -G, --max-gap-length=N       max gap size to fill in bp [800]; runtime increases\n"
"                               exponentially with respect to this parameter\n"
"  -j, --threads=N              use N parallel threads [1]\n"
"  -k, --kmer=N                 the size of a k-mer\n"
"  -b, --bloom-size=N           size of Bloom filter (e.g. '40G'). Required\n"
"                               when not using pre-built Bloom filter(s)\n"
"                               (-i option)\n"
"  -B, --max-branches=N         max branches in de Bruijn graph traversal;\n"
"                               use 'nolimit' for no limit [1000]\n"
"  -d, --dot-file=FILE          write graph traversals to a DOT file\n"
"  -e, --fix-errors             find and fix single-base errors when reads\n"
"                               have no kmers in bloom filter [disabled]\n"
"  -i, --input-bloom=FILE       load bloom filter from FILE\n"
"      --mask                   mask new and changed bases as lower case\n"
"      --no-mask                do not mask bases [default]\n"
"      --chastity               discard unchaste reads [default]\n"
"      --no-chastity            do not discard unchaste reads\n"
"      --trim-masked            trim masked bases from the ends of reads\n"
"      --no-trim-masked         do not trim masked bases from the ends\n"
"                               of reads [default]\n"
"  -m, --flank-mismatches=N     max mismatches between paths and flanks; use\n"
"                               'nolimit' for no limit [nolimit]\n"
"  -M, --max-mismatches=N       max mismatches between all alternate paths;\n"
"                               use 'nolimit' for no limit [nolimit]\n"
"  -n  --no-limits              disable all limits; equivalent to\n"
"                               '-B nolimit -m nolimit -M nolimit -P nolimit'\n"
"  -o, --output-prefix=FILE     prefix of output FASTA files [required]\n"
"  -P, --max-paths=N            merge at most N alternate paths; use 'nolimit'\n"
"                               for no limit [2]\n"
"  -q, --trim-quality=N         trim bases from the ends of reads whose\n"
"                               quality is less than the threshold\n"
"      --standard-quality       zero quality is `!' (33)\n"
"                               default for FASTQ and SAM files\n"
"      --illumina-quality       zero quality is `@' (64)\n"
"                               default for qseq and export files\n"
"  -r, --read-name=STR          only process reads with names that contain STR\n"
"  -s, --search-mem=N           mem limit for graph searches; multiply by the\n"
"                               number of threads (-j) to get the total mem used\n"
"                               for graph traversal [500M]\n"
"  -g, --gap-file=FILE          write sealed gaps to FILE\n"
"  -t, --trace-file=FILE        write graph search stats to FILE\n"
"  -v, --verbose                display verbose output\n"
"      --help                   display this help and exit\n"
"      --version                output version information and exit\n"
"\n"
" Deprecated Options:\n"
"\n"
"  -f, --min-frag=N             min fragment size in base pairs\n"
"  -F, --max-frag=N             max fragment size in base pairs\n"
"\n"
"   Note: --max-frag was formerly used to determine the maximum gap\n"
"   size that abyss-sealer would attempt to close, according to the formula\n"
"   max_gap_size = max_frag - 2 * flank_length, where flank_length is\n"
"   deteremined by the -L option.  --max-frag is kept only for backwards\n"
"   compatibility and is superceded by the more intuitive -G (--max-gap-length)\n"
"   option. Similarly, --min-frag determines the minimum gap size to close,\n"
"   according to the formula min_gap_size = min_frag - 2 * flank_length, where\n"
"   a negative gap size indicates an overlap between gap flanks.  Normally the\n"
"   user would not want to specify a minimum gap size and so it is recommended to\n"
"   leave --min-frag unset.\n"
"\n"
"Report bugs to <" PACKAGE_BUGREPORT ">.\n";

namespace opt {

	/** Length of flank. */
	int flankLength = 100;

	/** Max gap size to fill */
	unsigned maxGapLength = 800;

	/** scaffold file input. */
	static string inputScaffold;

	/** The number of parallel threads. */
	static unsigned threads = 1;

	/** The size of the bloom filter in bytes. */
	size_t bloomSize = 0;

	/** The maximum count value of the BLoom filter. */
	unsigned max_count = 2;

	/** Input read files are interleaved? */
	bool interleaved = false;

	/** Max active branches during de Bruijn graph traversal */
	unsigned maxBranches = 1000;

	/** multi-graph DOT file containing graph traversals */
	static string dotPath;

	/**
	 * Find and fix single base errors when a read has no
	 * kmers in the bloom filter.
	 */
	bool fixErrors = false;

	/** The size of a k-mer. */
	unsigned k;

	/** Vector of kmers. */
	vector<unsigned> kvector;

	/** Vector of Bloom filter paths corresponding to kmer values */
	vector<string> bloomFilterPaths;

	/** The minimum fragment size */
	unsigned minFrag = 0;

	/** The maximum fragment size */
	unsigned maxFrag = 0;

	/** Bloom filter input file */
	static string inputBloomPath;

	/** Max paths between left and right flanking sequences */
	unsigned maxPaths = 2;

	/** Prefix for output files */
	static string outputPrefix;

	/** Max mismatches allowed when building consensus seqs */
	unsigned maxMismatches = NO_LIMIT;

	/** Only process flanks that contain this substring. */
	static string readName;

	/** Max mem used per thread during graph traversal */
	static size_t searchMem = 500 * 1024 * 1024;

	/** Output file for graph search stats */
	static string tracefilePath;

	/** Output file for sealed gaps */
	static string gapfilePath;

	/** Mask bases not in flanks */
	static int mask = 0;

	/** Max mismatches between consensus and flanks */
	static unsigned maxFlankMismatches = NO_LIMIT;

	/** Output flanks files */
	static int printFlanks = 0;

	/** Output detailed stats */
	static int detailedStats = 0;
}

/** Counters */

struct Counters {
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
};

static const char shortopts[] = "S:L:b:B:d:ef:F:G:g:i:Ij:k:lm:M:no:P:q:r:s:t:v";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
	{ "detailed-stats",   no_argument, &opt::detailedStats, 1},
	{ "print-flanks",     no_argument, &opt::printFlanks, 1},
	{ "input-scaffold",   required_argument, NULL, 'S' },
	{ "flank-length",     required_argument, NULL, 'L' },
	{ "max-gap-length",   required_argument, NULL, 'G' },
	{ "bloom-size",       required_argument, NULL, 'b' },
	{ "max-branches",     required_argument, NULL, 'B' },
	{ "dot-file",         required_argument, NULL, 'd' },
	{ "fix-errors",       no_argument, NULL, 'e' },
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
	{ "flank-mismatches", required_argument, NULL, 'm' },
	{ "max-mismatches",   required_argument, NULL, 'M' },
	{ "max-paths",        required_argument, NULL, 'P' },
	{ "trim-quality",     required_argument, NULL, 'q' },
	{ "standard-quality", no_argument, &opt::qualityOffset, 33 },
	{ "illumina-quality", no_argument, &opt::qualityOffset, 64 },
	{ "read-name",        required_argument, NULL, 'r' },
	{ "search-mem",       required_argument, NULL, 's' },
	{ "trace-file",       required_argument, NULL, 't' },
	{ "gap-file",         required_argument, NULL, 'g' },
	{ "verbose",          no_argument, NULL, 'v' },
	{ "help",             no_argument, NULL, OPT_HELP },
	{ "version",          no_argument, NULL, OPT_VERSION },
	{ NULL, 0, NULL, 0 }
};

/** The coordinates of a feature. */
struct Coord
{
	int start;
	int end;

	Coord() { }
	Coord(int start, int end) : start(start), end(end) { }

	/** Return the size of this feature. */
	int size() const { return end - start; }
};

/** The coordinates of a gap and its flanking scaftigs. */
struct Gap
{
	/** The left flanking scaftig. */
	Coord left;

	/** The right flanking scaftig. */
	Coord right;

	Gap() { }

	Gap(int leftStart, int leftEnd, int rightStart, int rightEnd)
		: left(leftStart, leftEnd), right(rightStart, rightEnd) { }

	/** The start of the gap. */
	int gapStart() const { return left.end; }

	/** The end of the gap. */
	int gapEnd() const { return right.start; }

	/** The total size of the flanks plus the gap. */
	int totalSize() const { return right.end - left.start; }

	/** The size of the original gap. */
	int gapSize() const { return gapEnd() - gapStart(); }

	/** The size of the flanks. */
	int flankSize() const { return left.size() + right.size(); }
};

/** The coordinates of a gap, its flanking scaftigs, and the sequence that
 * fills it.
 */
struct ClosedGap : Gap
{
	/** The sequence that fills the gap. */
	std::string seq;

	ClosedGap() { }

	ClosedGap(Gap gap, std::string seq) : Gap(gap), seq(seq) { }
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

/**
 * Set the value for a commandline option, using "nolimit"
 * to represent NO_LIMIT.
 */
static inline void setMaxOption(unsigned& arg, istream& in)
{
	string str;
	getline(in, str);
	if (in && str == "nolimit") {
		arg = NO_LIMIT;
	} else {
		istringstream ss(str);
		ss >> arg;
		// copy state bits (fail, bad, eof) to
		// original stream
		in.clear(ss.rdstate());
	}
}

string sizetToString (size_t a)
{
	ostringstream temp;
	temp << a;
	return temp.str();
}

string IntToString (int a)
{
	ostringstream temp;
	temp<<a;
	return temp.str();
}

// returns merged sequence resulting from Konnector
template <typename Graph>
string merge(const Graph& g,
	unsigned k,
	const Gap& gap,
	const FastaRecord &read1,
	const FastaRecord &read2,
	const ConnectPairsParams& params,
	Counters& g_count,
	ofstream& traceStream)
{
	ConnectPairsResult result = connectPairs(k, read1, read2, g, params);
	ostringstream ss;
	ss << result.readNamePrefix << '_' << gap.gapStart() << '_' << gap.gapSize();
	result.readNamePrefix = ss.str();

	if (!opt::tracefilePath.empty())
#pragma omp critical(tracefile)
	{
			traceStream << result;
			assert_good(traceStream, opt::tracefilePath);
	}

	vector<FastaRecord>& paths = result.mergedSeqs;
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
#pragma omp atomic
					++g_count.tooManyReadMismatches;
			}
			else if (paths.size() > 1) {
#pragma omp atomic
				++g_count.multiplePaths;
			}
			else {
#pragma omp atomic
				++g_count.uniquePath;
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

	if (result.pathResult == FOUND_PATH) {
		if  (result.pathMismatches > params.maxPathMismatches)
			return "";
		else if (result.mergedSeqs.size() > 1)
			return result.consensusSeq;
		else
			return result.mergedSeqs.front();
	}
	else
		return "";
}

void printLog(ofstream &logStream, const string &output) {
#pragma omp critical(logStream)
	logStream << output;
	if (opt::verbose > 0)
#pragma omp critical(cerr)
		cerr << output;
}

void insertIntoScaffold(ofstream &scaffoldStream,
	ofstream &mergedStream,
	FastaRecord &record,
	map<string, map<int, ClosedGap> > &allmerged,
	unsigned &gapsclosedfinal)
{
	map<string, map<int, ClosedGap> >::iterator scaf_it;
	map<int, ClosedGap>::reverse_iterator pos_it;

	scaf_it = allmerged.find(record.id);
	scaffoldStream << ">" << record.id << " " << record.comment << endl;
	string modifiedSeq = record.seq;
	if (scaf_it != allmerged.end()) {
		for (pos_it = allmerged[record.id].rbegin();
			pos_it != allmerged[record.id].rend();
			pos_it++)
		{
			const ClosedGap& closedGap = pos_it->second;
			int newseqsize = pos_it->second.seq.length() - closedGap.flankSize();
			modifiedSeq.replace(
				closedGap.left.start,
				closedGap.totalSize(),
				closedGap.seq
			);
			mergedStream << ">" << record.id << "_" << pos_it->first << "_" << newseqsize <<  endl;
			mergedStream << pos_it->second.seq << endl;
			gapsclosedfinal++;
		}
		scaffoldStream << modifiedSeq << endl;
	}
	else {
		scaffoldStream << record.seq  << endl;
	}

}

std::pair<FastaRecord, FastaRecord>
makePseudoReads(std::string seq, FastaRecord record, const Gap& gap)
{
	std::pair<FastaRecord, FastaRecord> scaftigs;
	// extract left flank
	std::string leftflank = seq.substr(gap.left.start, gap.left.size());
	std::transform(leftflank.begin(), leftflank.end(), leftflank.begin(), ::toupper);
	scaftigs.first.seq = leftflank;
	scaftigs.first.id = record.id + "/1";

	// extract right flank
	scaftigs.second.id = record.id + "/2";
	std::string rightflank = seq.substr(gap.right.start, gap.right.size());
	std::transform(rightflank.begin(), rightflank.end(), rightflank.begin(), ::toupper);
	scaftigs.second.seq = reverseComplement(rightflank);

	return scaftigs;
}

template <typename Graph>
void kRun(const ConnectPairsParams& params,
	unsigned k,
	const Graph& g,
	map<string, map<int, ClosedGap> > &allmerged,
	map<FastaRecord, map<FastaRecord, Gap> > &flanks,
	unsigned &gapsclosed,
	ofstream &logStream,
	ofstream &traceStream,
	ofstream &gapStream)
{
	map<FastaRecord, map<FastaRecord, Gap> >::iterator read1_it;
	map<FastaRecord, Gap>::iterator read2_it;
	unsigned uniqueGapsClosed = 0;

	Counters g_count;
	g_count.noStartOrGoalKmer = 0;
	g_count.noPath = 0;
	g_count.uniquePath = 0;
	g_count.multiplePaths = 0;
	g_count.tooManyPaths = 0;
	g_count.tooManyBranches = 0;
	g_count.tooManyMismatches = 0;
	g_count.tooManyReadMismatches = 0;
	g_count.containsCycle = 0;
	g_count.exceededMemLimit = 0;
	g_count.traversalMemExceeded = 0;
	g_count.readPairsProcessed = 0;
	g_count.readPairsMerged = 0;
	g_count.skipped = 0;

	printLog(logStream, "Flanks inserted into k run = " + IntToString(flanks.size()) + "\n");

	int counter = 0;
	vector<map<FastaRecord, map<FastaRecord, Gap> >::iterator> flanks_closed;
#pragma omp parallel private(read1_it, read2_it) firstprivate(counter)
	for (read1_it = flanks.begin(); read1_it != flanks.end(); ++read1_it) {
		FastaRecord read1 = read1_it->first;
		bool success = false;
		for (read2_it = flanks[read1].begin(); read2_it != flanks[read1].end(); ++read2_it, ++counter) {
#if _OPENMP
			if (counter % omp_get_num_threads() != omp_get_thread_num())
				continue;
#endif
			FastaRecord read2 = read2_it->first;

			int startposition = read2_it->second.gapStart();
			string tempSeq = merge(g, k, read2_it->second, read1, read2, params, g_count, traceStream);
			if (!tempSeq.empty()) {
				success = true;
#pragma omp critical (allmerged)
				allmerged[read1.id.substr(0,read1.id.length()-2)][startposition]
					= ClosedGap(read2_it->second, tempSeq);
#pragma omp atomic
				++uniqueGapsClosed;
#pragma omp critical (gapsclosed)
				if (++gapsclosed % 100 == 0)
					printLog(logStream, IntToString(gapsclosed) + " gaps closed so far\n");

				if (!opt::gapfilePath.empty())
#pragma omp critical (gapStream)
					gapStream << ">" << read1.id.substr(0,read1.id.length()-2)
						  << "_" << read2_it->second.gapStart() << "-" << read2_it->second.gapEnd()
						  << " LN:i:" << tempSeq.length() << '\n'
						  << tempSeq << '\n';
			}
		}
		if (success) {
#pragma omp critical (flanks_closed)
			flanks_closed.push_back(read1_it);
		}
	}

	for (vector<map<FastaRecord, map<FastaRecord, Gap> >::iterator>::iterator it = flanks_closed.begin();
			it != flanks_closed.end(); ++it) {
		if (flanks.count((*it)->first) > 0)
		    flanks.erase(*it);
	}

	printLog(logStream, IntToString(uniqueGapsClosed) + " unique gaps closed for k" + IntToString(k) + "\n");

	printLog(logStream, "No start/goal kmer: "
		+ sizetToString(g_count.noStartOrGoalKmer) + "\n");
	printLog(logStream, "No path: "
		+ sizetToString(g_count.noPath) + "\n");
	printLog(logStream, "Unique path: "
		+ sizetToString(g_count.uniquePath) + "\n");
	printLog(logStream, "Multiple paths: "
		+ sizetToString(g_count.multiplePaths) + "\n");
	printLog(logStream, "Too many paths: "
		+ sizetToString(g_count.tooManyPaths) + "\n");
	printLog(logStream, "Too many branches: "
		+ sizetToString(g_count.tooManyBranches) + "\n");
	printLog(logStream, "Too many path/path mismatches: "
		+ sizetToString(g_count.tooManyMismatches) + "\n");
	printLog(logStream, "Too many path/read mismatches: "
		+ sizetToString(g_count.tooManyReadMismatches) + "\n");
	printLog(logStream, "Contains cycle: "
		+ sizetToString(g_count.containsCycle) + "\n");
	printLog(logStream, "Exceeded mem limit: "
		+ sizetToString(g_count.exceededMemLimit) + "\n");
	printLog(logStream, "Skipped: "
		+ sizetToString(g_count.skipped) + "\n");

	printLog(logStream, IntToString(flanks.size()) + " flanks left\n");

}

bool operator<(const FastaRecord& a, const FastaRecord& b)
{
	if (a.id == b.id)
		return a.seq < b.seq;
	else
		return a.id < b.id;
}

void findFlanks(FastaRecord &record,
	int flanklength,
	unsigned &gapnumber,
	map<FastaRecord, map<FastaRecord, Gap> > &flanks)
{
	const string& seq = record.seq;

	// Iterate over the gaps.
	for (size_t offset = 0;
			seq.string::find_first_of("Nn", offset) != string::npos;) {
#pragma omp atomic
		gapnumber++;
		size_t startposition = seq.string::find_first_of("Nn", offset);

		size_t endposition = seq.string::find_first_not_of("Nn", startposition);
		if (endposition == string::npos) {
			std::cerr << PROGRAM ": Warning: sequence ends with an N: " << record.id << "\n";
			break;
		}

		size_t right_end = seq.string::find_first_of("Nn", endposition);
		if (right_end == string::npos)
			right_end = seq.length();

		Gap gap(
			max((ssize_t)offset, (ssize_t)startposition - (ssize_t)flanklength),
			startposition,
			endposition,
			min(right_end, endposition + flanklength));

		std::pair<FastaRecord, FastaRecord> scaftigs =
			makePseudoReads(seq, record, gap);
#pragma omp critical (flanks)
		flanks[scaftigs.first][scaftigs.second] = gap;

		offset = endposition;
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
		  case 'S':
			arg >> opt::inputScaffold; break;
		  case 'L':
			arg >> opt::flankLength; break;
		  case 'G':
			arg >> opt::maxGapLength; break;
		  case 'b':
			opt::bloomSize = SIToBytes(arg); break;
		  case 'B':
			setMaxOption(opt::maxBranches, arg); break;
		  case 'd':
			arg >> opt::dotPath; break;
		  case 'e':
			opt::fixErrors = true; break;
		  case 'f':
			arg >> opt::minFrag; break;
		  case 'F':
			arg >> opt::maxFrag; break;
		  case 'i': {
			string tempPath;
			arg >> tempPath;
			opt::bloomFilterPaths.push_back(tempPath);
			opt::inputBloomPath = tempPath;
			break;
			}
		  case 'I':
			opt::interleaved = true; break;
		  case 'j':
			arg >> opt::threads; break;
		  case 'k': {
			unsigned tempK;
			arg >> tempK;
			opt::kvector.push_back(tempK);
			opt::k = tempK;
			break;
			}
		  case 'm':
			setMaxOption(opt::maxFlankMismatches, arg); break;
		  case 'n':
			opt::maxBranches = NO_LIMIT;
			opt::maxFlankMismatches = NO_LIMIT;
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
		  case 'g':
		    arg >> opt::gapfilePath; break;
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

	/* translate --max-frag to --max-gap-length for backwards compatibility */
	if (opt::maxFrag > 0) {
		if ((int)opt::maxFrag < 2 * opt::flankLength)
			opt::maxGapLength = 0;
		else
			opt::maxGapLength = opt::maxFrag - 2 * opt::flankLength;
	}

	if (opt::inputScaffold.empty()) {
		cerr << PROGRAM ": missing mandatory option `-S'\n";
		die = true;
	}

	if (opt::k == 0) {
		cerr << PROGRAM ": missing mandatory option `-k'\n";
		die = true;
	}

	if (opt::bloomFilterPaths.size() < opt::kvector.size()
		&& opt::bloomSize == 0)
	{
		cerr << PROGRAM ": missing mandatory option `-b' (Bloom filter size)\n"
			<< "Here are some guidelines for sizing the Bloom filter:\n"
			<< "  * E. coli (~5 Mbp genome), 615X coverage: -b500M\n"
			<< "  * S. cerevisiae (~12 Mbp genome), 25X coverage: -b500M\n"
			<< "  * C. elegans (~100 Mbp genome), 89X coverage: -b1200M\n"
			<< "  * H. sapiens (~3 Gbp genome), 71X coverage: -b40G\n";
		die = true;
	}

	if (opt::outputPrefix.empty()) {
		cerr << PROGRAM ": missing mandatory option `-o'\n";
		die = true;
	}

	if (opt::bloomFilterPaths.size() > opt::kvector.size()) {
		cerr << PROGRAM ": you must specify a k-mer size (-k) for each Bloom "
			" filter file (-i)\n";
		die = true;
	} else if (opt::bloomFilterPaths.size() < opt::kvector.size()
		&& argc - optind < 1) {
		cerr << PROGRAM ": missing input file arguments\n";
		die = true;
	} else if (opt::bloomFilterPaths.size() == opt::kvector.size()
		&& argc - optind > 0) {
		cerr << PROGRAM ": input FASTA/FASTQ args should be omitted when using "
			"pre-built Bloom filters (-i) for all k-mer sizes\n";
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

	ofstream gapStream;
	if (!opt::gapfilePath.empty()) {
		gapStream.open(opt::gapfilePath.c_str());
		assert(gapStream.is_open());
		assert_good(gapStream, opt::gapfilePath);
	}

	string logOutputPath(opt::outputPrefix);
	logOutputPath.append("_log.txt");
	ofstream logStream(logOutputPath.c_str());
	assert_good(logStream, logOutputPath);

	string scaffoldOutputPath(opt::outputPrefix);
	scaffoldOutputPath.append("_scaffold.fa");
	ofstream scaffoldStream(scaffoldOutputPath.c_str());
	assert_good(scaffoldStream, scaffoldOutputPath);

	string mergedOutputPath(opt::outputPrefix);
	mergedOutputPath.append("_merged.fa");
	ofstream mergedStream(mergedOutputPath.c_str());
	assert_good(mergedStream, mergedOutputPath);

	printLog(logStream, "Finding flanks\n");

	ConnectPairsParams params;

	params.minMergedSeqLen = opt::minFrag;
	params.maxMergedSeqLen = opt::maxGapLength + 2 * opt::flankLength;
	params.maxPaths = opt::maxPaths;
	params.maxBranches = opt::maxBranches;
	params.maxPathMismatches = opt::maxMismatches;
	params.maxReadMismatches = opt::maxFlankMismatches;
	/*
	 * search from the end of the scaffold flank until
	 * we find this many matches in a row.
	 */
	params.kmerMatchesThreshold = 1;
	params.fixErrors = opt::fixErrors;
	params.maskBases = opt::mask;
	params.memLimit = opt::searchMem;
	params.dotPath = opt::dotPath;
	params.dotStream = opt::dotPath.empty() ? NULL : &dotStream;

	map<FastaRecord, map<FastaRecord, Gap> > flanks;
	const char* scaffoldInputPath = opt::inputScaffold.c_str();
	FastaReader reader1(scaffoldInputPath, FastaReader::FOLD_CASE);
	unsigned gapsfound = 0;
	string temp;

#pragma omp parallel
	for (FastaRecord record;;) {
		bool good;
#pragma omp critical(reader1)
		good = reader1 >> record;
		if (good) {
			findFlanks(record, opt::flankLength, gapsfound, flanks);
		}
		else {
			break;
		}
	}

	temp = IntToString(gapsfound) + " gaps found\n";
	printLog(logStream, temp);
	temp = IntToString((int)flanks.size()) + " flanks extracted\n\n";
	printLog(logStream, temp);

	if (opt::printFlanks > 0) {
		map<FastaRecord, map<FastaRecord, Gap> >::iterator read1_it;
		map<FastaRecord, Gap>::iterator read2_it;

		string read1OutputPath(opt::outputPrefix);
		read1OutputPath.append("_flanks_1.fa");
		ofstream read1Stream(read1OutputPath.c_str());
		assert_good(read1Stream, read1OutputPath);

		string read2OutputPath(opt::outputPrefix);
		read2OutputPath.append("_flanks_2.fa");
		ofstream read2Stream(read2OutputPath.c_str());
		assert_good(read2Stream, read2OutputPath);

		for (read1_it = flanks.begin(); read1_it != flanks.end(); read1_it++) {
			FastaRecord read1 = read1_it->first;
			for (read2_it = flanks[read1].begin(); read2_it != flanks[read1].end(); read2_it++) {
				FastaRecord read2 = read2_it->first;
				const Gap& gap = read2_it->second;

				read1Stream << ">" << read1.id.substr(0,read2.id.length()-2)
					<< "_" << gap.gapStart() << "_" << gap.gapSize() << "/1\n";
				read1Stream << read1.seq << endl;

				read2Stream << ">" << read2.id.substr(0,read2.id.length()-2)
					<< "_" << gap.gapStart() << "_" << gap.gapSize() << "/2\n";
				read2Stream << read2.seq << endl;
			}
		}
		assert_good(read1Stream, read1OutputPath.c_str());
		read1Stream.close();
		assert_good(read2Stream, read2OutputPath.c_str());
		read2Stream.close();
	}

	/** map for merged sequence resutls */
	map<string, map<int, ClosedGap> > allmerged;
	unsigned gapsclosed=0;

	for (unsigned i = 0; i<opt::kvector.size(); i++) {
		opt::k = opt::kvector.at(i);
		Kmer::setLength(opt::k);

		BloomFilter* bloom;
		CascadingBloomFilter* cascadingBloom = NULL;

		if (!opt::bloomFilterPaths.empty() && i < opt::bloomFilterPaths.size()) {

			temp = "Loading bloom filter from `" + opt::bloomFilterPaths.at(i) + "'...\n";
			printLog(logStream, temp);

			bloom = new BloomFilter();

			const char* inputPath = opt::bloomFilterPaths.at(i).c_str();
			ifstream inputBloom(inputPath, ios_base::in | ios_base::binary);
			assert_good(inputBloom, inputPath);
			inputBloom >> *bloom;
			assert_good(inputBloom, inputPath);
			inputBloom.close();
		} else {
			printLog(logStream, "Building bloom filter\n");

			size_t bits = opt::bloomSize * 8 / 2;
			cascadingBloom = new CascadingBloomFilter(bits, opt::max_count);
#ifdef _OPENMP
			ConcurrentBloomFilter<CascadingBloomFilter> cbf(*cascadingBloom, 1000);
			for (int i = optind; i < argc; i++)
				Bloom::loadFile(cbf, opt::k, argv[i], 0 /*opt::verbose*/);
#else
			for (int i = optind; i < argc; i++)
				Bloom::loadFile(*cascadingBloom, opt::k, argv[i], 0 /* opt::verbose*/);
#endif
			bloom = &cascadingBloom->getBloomFilter(opt::max_count - 1);
		}

		assert(bloom != NULL);

		if (opt::verbose)
			cerr << "Bloom filter FPR: " << setprecision(3)
				<< 100 * bloom->FPR() << "%\n";

		DBGBloom<BloomFilter> g(*bloom);

		temp = "Starting K run with k = " + IntToString(opt::k) + "\n";
		printLog(logStream, temp);

		kRun(params, opt::k, g, allmerged, flanks, gapsclosed, logStream, traceStream, gapStream);

		temp = "k" + IntToString(opt::k) + " run complete\n"
				+ "Total gaps closed so far = " + IntToString(gapsclosed) + "\n\n";
		printLog(logStream, temp);

		if (!opt::inputBloomPath.empty())
			delete bloom;
		else
			delete cascadingBloom;
	}

	printLog(logStream, "K sweep complete\nCreating new scaffold with gaps closed...\n");

	map<string, map<int, map<string, string> > >::iterator scaf_it;
	map<int, map<string, string> >::reverse_iterator pos_it;
	FastaReader reader2(scaffoldInputPath, FastaReader::FOLD_CASE);
	unsigned gapsclosedfinal = 0;

	/** creating new scaffold with gaps closed */
	for (FastaRecord record;;) {
		bool good;
		good = reader2 >> record;
		if (good) {
			insertIntoScaffold(scaffoldStream, mergedStream, record, allmerged, gapsclosedfinal);
		}
		else
			break;
	}
	printLog(logStream, "New scaffold complete\n");
	printLog(logStream, "Gaps closed = " + IntToString(gapsclosed) + "\n");
	logStream << (float)100 * gapsclosed / gapsfound << "%\n\n";
	if (opt::verbose > 0) {
		cerr << (float)100 * gapsclosed / gapsfound << "%\n\n";
	}

	assert_good(scaffoldStream, scaffoldOutputPath.c_str());
	scaffoldStream.close();
	assert_good(mergedStream, mergedOutputPath.c_str());
	mergedStream.close();
	assert_good(logStream, logOutputPath.c_str());
	logStream.close();

	if (!opt::dotPath.empty()) {
		assert_good(dotStream, opt::dotPath);
		dotStream.close();
	}

	if (!opt::tracefilePath.empty()) {
		assert_good(traceStream, opt::tracefilePath);
		traceStream.close();
	}

	if (!opt::gapfilePath.empty()) {
		assert_good(gapStream, opt::gapfilePath);
		gapStream.close();
	}

	return 0;
}
