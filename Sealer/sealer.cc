/**
 * Connect pairs using a Bloom filter de Bruijn graph
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

//#include "Uncompress.h"

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
"Usage: " PROGRAM " -k <kmer size> -k <kmer size>... -o <output_prefix> -S <path to scaffold file> [options]... <reads1> [reads2]...\n"
"i.e. abyss-sealer -k30 -k40 -k60 -k65 -o test -S scaffold.fa read1.fa read2.fa\n\n"
"Close gaps by using left and right flanks as reads for Konnector\n"
"and performing a K sweep of supplied K values..\n"
"\n"
" Options:\n"
"\n"
"      --detailed-stats		outputs detailed stats\n"
"      --print-flanks		outputs flank files\n"
"      --no-cascade		do not use cascading bloom filter [default]\n"
"  -S, --input-scaffold=FILE  	load scaffold from FILE\n"
"      --cascade			build cascading bloom filter\n"
"  -L, --flank-length=N	      	length of flanks to be used as pseudoreads [100]\n"
"  -D, --flank-distance=N     	distance of flank from gap [0]\n"
"  -j, --threads=N            	use N parallel threads [1]\n"
"  -k, --kmer=N			the size of a k-mer\n"
"  -b, --bloom-size=N        	size of bloom filter [500M]\n"
"  -B, --max-branches=N    	max branches in de Bruijn graph traversal;\n"
"                             	use 'nolimit' for no limit [350]\n"
"  -d, --dot-file=FILE        	write graph traversals to a DOT file\n"
"  -e, --fix-errors           	find and fix single-base errors when reads\n"
"                             	have no kmers in bloom filter [disabled]\n"
"  -f, --min-frag=N           	min fragment size in base pairs [0]\n"
"  -F, --max-frag=N           	max fragment size in base pairs [1000]\n"
"  -i, --input-bloom=FILE     	load bloom filter from FILE\n"
"  -I, --interleaved          	input reads files are interleaved\n"
"      --mask                 	mask new and changed bases as lower case\n"
"      --no-mask              	do not mask bases [default]\n"
"      --chastity             	discard unchaste reads [default]\n"
"      --no-chastity          	do not discard unchaste reads\n"
"      --trim-masked          	trim masked bases from the ends of reads\n"
"      --no-trim-masked       	do not trim masked bases from the ends\n"
"                             	of reads [default]\n"
"  -l, --long-search          	start path search as close as possible\n"
"                             	to the beginnings of reads. Takes more time\n"
"                             	but improves results when bloom filter false\n"
"                             	positive rate is high [disabled]\n"
"  -m, --read-mismatches=N    	max mismatches between paths and reads; use\n"
"                             	'nolimit' for no limit [nolimit]\n"
"  -M, --max-mismatches=N     	max mismatches between all alternate paths;\n"
"                             	use 'nolimit' for no limit [2]\n"
"  -n  --no-limits            	disable all limits; equivalent to\n"
"                             	'-B nolimit -m nolimit -M nolimit -P nolimit'\n"
"  -o, --output-prefix=FILE   	prefix of output FASTA files [required]\n"
"  -P, --max-paths=N          	merge at most N alternate paths; use 'nolimit'\n"
"                             	for no limit [2]\n"
"  -q, --trim-quality=N       	trim bases from the ends of reads whose\n"
"                             	quality is less than the threshold\n"
"      --standard-quality     	zero quality is `!' (33)\n"
"                             	default for FASTQ and SAM files\n"
"      --illumina-quality     	zero quality is `@' (64)\n"
"                             	default for qseq and export files\n"
"  -r, --read-name=STR        	only process reads with names that contain STR\n"
"  -s, --search-mem=N         	mem limit for graph searches; multiply by the\n"
"                             	number of threads (-j) to get the total mem used\n"
"                             	for graph traversal [500M]\n"
"  -t, --trace-file=FILE      	write graph search stats to FILE\n"
"  -v, --verbose              	display verbose output\n"
"      --help                 	display this help and exit\n"
"      --version              	output version information and exit\n"
"\n"
"Report bugs to <" PACKAGE_BUGREPORT ">.\n";

const unsigned g_progressStep = 1000;

namespace opt {

	/** Length of flank. */
	int flankLength = 100;

	/** Distance of flank from gap. */
	unsigned flankDistance = 0;

	/** scaffold file input. */
	static string inputScaffold;

	/** The number of parallel threads. */
	static unsigned threads = 1;

	/** The size of the bloom filter in bytes. */
	size_t bloomSize = 500 * 1024 * 1024;

	/** Input read files are interleaved? */
	bool interleaved = false;

	/**
	 * Choose start/goal kmers for path search as close as
	 * possible to beginning (5' end) of reads. Improves
	 * results when bloom filter FPR is high.
	 */
	bool longSearch = false;

	/** Max active branches during de Bruijn graph traversal */
	unsigned maxBranches = 350;

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

	/** Use cascading bloom filter in addition to plain bloom filter */
	static int cascade = 0;

	/** Output flanks files */
	static int printFlanks = 0;

	/** Output detailed stats */
	static int detailedStats = 0;
}

/** Counters */

//static struct {
//	size_t noStartOrGoalKmer;
//	size_t noPath;
//	size_t uniquePath;
//	size_t multiplePaths;
//	size_t tooManyPaths;
//	size_t tooManyBranches;
//	size_t tooManyMismatches;
//	size_t tooManyReadMismatches;
//	size_t containsCycle;
//	size_t exceededMemLimit;
//	size_t traversalMemExceeded;
//	size_t readPairsProcessed;
//	size_t readPairsMerged;
//	size_t skipped;
//} g_count;

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

//static Counters g_count;


static const char shortopts[] = "S:L:D:b:B:d:ef:F:i:Ij:k:lm:M:no:P:q:r:s:t:v";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
	{ "detailed-stats",   no_argument, &opt::detailedStats, 1},
	{ "cascade",          no_argument, &opt::cascade, 1 },
	{ "print-flanks",     no_argument, &opt::printFlanks, 1},
	{ "no-cascade",       no_argument, &opt::cascade, 0 },
	{ "input-scaffold",   required_argument, NULL, 'S' },
	{ "flank-length",     required_argument, NULL, 'L' },
	{ "flank-distance",   required_argument, NULL, 'D' },
	{ "bloom-size",       required_argument, NULL, 'b' },
	{ "max-branches",     required_argument, NULL, 'B' },
	{ "dot-file",         required_argument, NULL, 'd' },
	{ "fix-errors",       no_argument, NULL, 'e' },
	{ "min-frag",         required_argument, NULL, 'f' },
	{ "max-frag",         required_argument, NULL, 'F' },
	{ "input-bloom",      required_argument, NULL, 'i' },
	{ "interleaved",      no_argument, NULL, 'I' },
	{ "long-search",      no_argument, NULL, 'l' },
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



string reversecompliment(string str) {
        map<string, string> dict;
        dict["A"] = "T";
        dict["T"] = "A";
        dict["G"] = "C";
        dict["C"] = "G";
        dict["M"] = "K";
        dict["R"] = "Y";
        dict["W"] = "W";
        dict["S"] = "S";
        dict["Y"] = "R";
        dict["V"] = "B";
        dict["K"] = "M";
        dict["H"] = "D";
        dict["D"] = "H";
        dict["B"] = "V";
        dict["N"] = "N";

        unsigned base = 0;
        string result = "";
        while (base <= str.length()) {
                result = result + dict[str.substr(base, 1)];
                base++;
        }
        reverse(result.begin(), result.end());
        return result;
}

string IntToString (int a)
{
    ostringstream temp;
    temp<<a;
    return temp.str();
}

int StringToInt(string num) {
	stringstream str(num);
	int x;
	str >> x;
	return x;
}

// returns merged sequence resulting from Konnector
template <typename Graph>
string merge(const Graph& g,
	unsigned k,
	FastaRecord read1,
	FastaRecord read2,
	const ConnectPairsParams& params,
	Counters& g_count)
{
	ConnectPairsResult result = connectPairs(k, read1, read2, g, params);

	vector<FastaRecord>& paths = result.mergedSeqs;
	switch (result.pathResult) {

		case NO_PATH:
			assert(paths.empty());
			if (result.foundStartKmer && result.foundGoalKmer)
				++g_count.noPath;
			else {
				++g_count.noStartOrGoalKmer;
			}
			break;

		case FOUND_PATH:
			assert(!paths.empty());
			if (result.pathMismatches > params.maxPathMismatches ||
					result.readMismatches > params.maxReadMismatches) {
				if (result.pathMismatches > params.maxPathMismatches)
					++g_count.tooManyMismatches;
				else
					++g_count.tooManyReadMismatches;
			}
			else if (paths.size() > 1) {
				++g_count.multiplePaths;
			}
			else {
				++g_count.uniquePath;
			}
			break;

		case TOO_MANY_PATHS:
			++g_count.tooManyPaths;
			break;

		case TOO_MANY_BRANCHES:
			++g_count.tooManyBranches;
			break;

		case PATH_CONTAINS_CYCLE:
			++g_count.containsCycle;
			break;

		case EXCEEDED_MEM_LIMIT:
			++g_count.exceededMemLimit;
			break;
	}

	if (result.pathResult == FOUND_PATH) {
		if (result.mergedSeqs.size() > 1)
			return result.consensusSeq;
		else
			return result.mergedSeqs.front();
	}
	else
		return "";


	vector<FastaRecord>& paths = result.mergedSeqs;
	switch (result.pathResult) {

		case NO_PATH:
			assert(paths.empty());
			if (result.foundStartKmer && result.foundGoalKmer)
				++g_count.noPath;
			else {
				++g_count.noStartOrGoalKmer;
			}
			break;

		case FOUND_PATH:
			assert(!paths.empty());
			if (result.pathMismatches > params.maxPathMismatches ||
					result.readMismatches > params.maxReadMismatches) {
				if (result.pathMismatches > params.maxPathMismatches)
					++g_count.tooManyMismatches;
				else
					++g_count.tooManyReadMismatches;
			}
			else if (paths.size() > 1) {
				++g_count.multiplePaths;
			}
			else {
				++g_count.uniquePath;
			}
			break;

		case TOO_MANY_PATHS:
			++g_count.tooManyPaths;
			break;

		case TOO_MANY_BRANCHES:
			++g_count.tooManyBranches;
			break;

		case PATH_CONTAINS_CYCLE:
			++g_count.containsCycle;
			break;

		case EXCEEDED_MEM_LIMIT:
			++g_count.exceededMemLimit;
			break;
	}
}

void insertIntoScaffold(ofstream &scaffoldStream,
	ofstream &mergedStream,
	FastaRecord record,
	map<string, map<int, map<string, string> > > &allmerged,
	unsigned &gapsclosedfinal)
{
	map<string, map<int, map<string, string> > >::iterator scaf_it;
	map<int, map<string, string> >::reverse_iterator pos_it;

	scaf_it = allmerged.find(record.id);
	scaffoldStream << ">" << record.id << endl;
	string modifiedSeq = record.seq;
	if (scaf_it != allmerged.end()) {
		for (pos_it = allmerged[record.id].rbegin();
			pos_it != allmerged[record.id].rend();
			pos_it++)
		{
			modifiedSeq.replace(
				pos_it->first - opt::flankLength,
				StringToInt(pos_it->second["gap"]) + (opt::flankLength * 2),
				pos_it->second["seq"]
			);
			mergedStream << ">" << record.id << " gap closed"  << endl;
                        mergedStream << pos_it->second["seq"] << endl;
			gapsclosedfinal++;
		}
		scaffoldStream << modifiedSeq << endl;
	}
	else {
		scaffoldStream << record.seq  << endl;
	}

}

void makePseudoReads(
	FastaRecord &read1,
	FastaRecord &read2,
	int startposition,
	int endposition,
	string seq,
	int flanklength,
	FastaRecord record)
{
	// extract left flank
	string leftflank = seq.substr(startposition - flanklength, flanklength);
	transform(leftflank.begin(), leftflank.end(), leftflank.begin(), ::toupper); // requires <algorithm>
	read1.seq = leftflank;
	read1.id = record.id + "/1";

	// extract right flank
	read2.id = record.id + "/2";
	string rightflank = seq.substr(endposition, flanklength);
	transform(rightflank.begin(), rightflank.end(), rightflank.begin(), ::toupper);
	read2.seq = reversecompliment(rightflank);
}

template <typename Graph>
void kRun(const ConnectPairsParams& params,
	unsigned k,
	const Graph& g,
	map<string, map<int, map<string, string> > > &allmerged,
	map<FastaRecord, map<FastaRecord, map<string, int> > > &flanks,
	unsigned &gapsclosed)
{
	map<FastaRecord, map<FastaRecord, map<string, int> > >:: iterator read1_it;
	map<FastaRecord, map<string, int> >::iterator read2_it;
	unsigned uniqueGapsClosed = 0;
	bool success;
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


	if (opt::verbose > 0)
		cerr << "Reads inserted into k run = " << flanks.size() << endl;

	for (read1_it = flanks.begin(); read1_it != flanks.end();) {
		success = false;
		FastaRecord read1 = read1_it->first;
		for (read2_it = flanks[read1].begin(); read2_it != flanks[read1].end(); read2_it++) {
			FastaRecord read2 = read2_it->first;

			int startposition = read2_it->second["startposition"];
			int endposition = read2_it->second["endposition"];
			string tempSeq;
			tempSeq = merge(g, k, read1, read2, params, g_count);

			if (!tempSeq.empty()) {
				success = true;
				allmerged[read1.id.substr(0,read1.id.length()-2)][startposition]["gap"]
					= IntToString(endposition - startposition);
				allmerged[read1.id.substr(0,read1.id.length()-2)][startposition]["seq"]
					= tempSeq;
				gapsclosed++;
				uniqueGapsClosed++;
				if (gapsclosed % 100 == 0 && opt::verbose > 0)
					cerr << gapsclosed << " gaps closed so far" << endl;
			}
		}
		if (success)
			flanks.erase(read1_it++);
		else
			read1_it++;
	}

	if (opt::verbose > 0) {
		cerr << uniqueGapsClosed << " unique gaps closed for k"
			<< k << endl;

		if (opt::detailedStats > 0) {
			cerr <<
			//	"Processed " << initialFlankNumber << " pseudoreads\n"
			//	"Merged (Unique path + Multiple paths): "
			//		<< g_count.uniquePath + g_count.multiplePaths
			//		<< /*" (" << setprecision(3) <<  (float)100
			//		    * (g_count.uniquePath + g_count.multiplePaths) /
			//		   initialFlankNumber
			//		<< "%)*/"\n"
				"No start/goal kmer: " << g_count.noStartOrGoalKmer
					<< /*" (" << setprecision(3) << (float)100
						* g_count.noStartOrGoalKmer / initialFlankNumber
					<< "%)*/"\n"
				"No path: " << g_count.noPath
					<< /*" (" << setprecision(3) << (float)100
						* g_count.noPath / initialFlankNumber
					<< "%)*/"\n"
				"Unique path: " << g_count.uniquePath
					<</* " (" << setprecision(3) << (float)100
						* g_count.uniquePath / initialFlankNumber
					<< "%)*/"\n"
				"Multiple paths: " << g_count.multiplePaths
					<</* " (" << setprecision(3) << (float)100
						* g_count.multiplePaths / initialFlankNumber
					<< "%)*/"\n"
				"Too many paths: " << g_count.tooManyPaths
					<</* " (" << setprecision(3) << (float)100
						* g_count.tooManyPaths / initialFlankNumber
					<< "%)*/"\n"
				"Too many branches: " << g_count.tooManyBranches
					<</* " (" << setprecision(3) << (float)100
						* g_count.tooManyBranches / initialFlankNumber
					<< "%)*/"\n"
				"Too many path/path mismatches: " << g_count.tooManyMismatches
					<</* " (" << setprecision(3) << (float)100
						* g_count.tooManyMismatches / initialFlankNumber
					<< "%)*/"\n"
				"Too many path/read mismatches: " << g_count.tooManyReadMismatches
					<< /*" (" << setprecision(3) << (float)100
						* g_count.tooManyReadMismatches / initialFlankNumber
					<< "%)*/"\n"
				"Contains cycle: " << g_count.containsCycle
					<< /*" (" << setprecision(3) << (float)100
						* g_count.containsCycle / initialFlankNumber
					<< "%)*/"\n"
				"Exceeded mem limit: " << g_count.exceededMemLimit
					<</* " (" << setprecision(3) << (float)100
						* g_count.exceededMemLimit / initialFlankNumber
					<< "%)*/"\n"
				"Skipped: " << g_count.skipped
					<</* " (" << setprecision(3) << (float)100
						* g_count.skipped / initialFlankNumber
					<< "%)*/"\n"
			/*	"Bloom filter FPR: " << setprecision(3) << 100 * bloom->FPR()
					<< "%\n"*/;
		}
		cerr << flanks.size() << " reads left" << endl;
	}

}

typedef map<string, int> property_map;
typedef map<FastaRecord, property_map> read2_map;
typedef map<FastaRecord, read2_map> read1_map;

bool operator<(const FastaRecord& a, const FastaRecord& b)
{
	if (a.id == b.id)
		return a.seq < b.seq;
	else
		return a.id < b.id;
}

void findFlanks(FastaRecord record,
	int flanklength,
	unsigned &gapnumber,
	read1_map &flanks)
{
        int offset = 0;
	int endposition = 0;
       	int startposition = 0;
        string seq = record.seq;
	unsigned seqsize = seq.length();


	// finds next gap. If no more gaps, while loop ends.
	while (seq.string::find_first_of("Nn", offset) != string::npos) {
		gapnumber++;
		startposition  = seq.string::find_first_of("Nn", offset);

		// if gap is at the end of the sequence (no right flank), the while loop ends.
		if (seq.string::find_first_not_of("Nn", startposition) == string::npos)
			break;
		else {
			// position of the first non-N character
			endposition = seq.string::find_first_not_of("Nn", startposition);
		}
		// if flanks are not long enough, move to next gap.
		if (startposition-flanklength<0 || unsigned(endposition + flanklength) > seqsize) {
			offset = endposition;
			continue;
		}
		// if flanks contain at least one N, move to next gap.
		else if (seq.substr(startposition - flanklength, flanklength).string::find_first_of("Nn")
				!= string::npos
			or seq.substr(endposition, flanklength).string::find_first_of("Nn")
				!= string::npos)
		{
			offset = endposition;
			continue;
		}

		FastaRecord read1, read2;
		makePseudoReads(read1, read2, startposition, endposition, seq, flanklength, record);

		flanks[read1][read2]["startposition"] = startposition;
		flanks[read1][read2]["endposition"] = endposition;

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
		  case 'D':
			arg >> opt::flankDistance; break;
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
		  case 'l':
			opt::longSearch = true; break;
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

	if (opt::inputScaffold.empty()) {
		cerr << PROGRAM ": missing mandatory option `-S'\n";
		die = true;
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
	/*
	BloomFilterBase* bloom = NULL;

	if (!opt::inputBloomPath.empty()) {

		if (opt::verbose)
			std::cerr << "Loading bloom filter from `"
				<< opt::inputBloomPath << "'...\n";

		const char* inputPath = opt::inputBloomPath.c_str();
		ifstream inputBloom(inputPath, ios_base::in | ios_base::binary);
		assert_good(inputBloom, inputPath);
		inputBloom >> bloom;
		assert_good(inputBloom, inputPath);
		inputBloom.close();

	} else {
		// Specify bloom filter size in bits.
		size_t bits = opt::bloomSize * 8 / 2;
		CascadingBloomFilter tempBloom(bits);
#ifdef_OPENMP
		ConcurrentBloomFilter<CascadingBloomfilter> cbf(tempBloom, 1000);
		for(int i = optind; i < argc; i++)
			Bloom::loadFile(cbf, opt::k, string(argv[i]), opt::verbose);
#else
		for (int i = optind; i < argc; i++)
			Bloom::loadFile(tempBloom, opt::k, string(argv[i]), opt::verbose);
#endif
		bloom = tempBloom.getBloomfilter(tempBloom.MAX_COUNT-1);
	}

	//cout << *bloom;
	//exit(EXIT_SUCCESS);

	if (opt::verbose)
		cerr << "Bloom filter FPR: " << setprecision(3)
			<< 100 * bloom->FPR() << "%\n";
*/

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

//	DBGBloom<BloomFilter> g(bloom);

	string scaffoldOutputPath(opt::outputPrefix);
	scaffoldOutputPath.append("_scaffold.fa");
	ofstream scaffoldStream(scaffoldOutputPath.c_str());
	assert_good(scaffoldStream, scaffoldOutputPath);

	string mergedOutputPath(opt::outputPrefix);
	mergedOutputPath.append("_merged.fa");
	ofstream mergedStream(mergedOutputPath.c_str());
	assert_good(mergedStream, mergedOutputPath);

	if (opt::verbose > 0)
		cerr << "Finding flanks\n";

	ConnectPairsParams params;

	params.minMergedSeqLen = opt::minFrag;
	params.maxMergedSeqLen = opt::maxFrag;
	params.maxPaths = opt::maxPaths;
	params.maxBranches = opt::maxBranches;
	params.maxPathMismatches = opt::maxMismatches;
	params.maxReadMismatches = opt::maxReadMismatches;
	params.fixErrors = opt::fixErrors;
	params.longSearch = opt::longSearch;
	params.maskBases = opt::mask;
	params.memLimit = opt::searchMem;
	params.dotPath = opt::dotPath;
	params.dotStream = opt::dotPath.empty() ? NULL : &dotStream;

	map<FastaRecord, map<FastaRecord, map<string, int> > > flanks;
	const char* scaffoldInputPath = opt::inputScaffold.c_str();
	FastaReader reader1(scaffoldInputPath, FastaReader::FOLD_CASE);
	unsigned gapsfound = 0;

	for (FastaRecord record;;) {
               	bool good;
               	good = reader1 >> record;
               	if (good) {
			findFlanks(record, opt::flankLength, gapsfound, flanks);
		}
               	else {
			if (opt::verbose > 0) {
				cerr << gapsfound << " gaps found" << endl;
				cerr << flanks.size() << " flanks extracted\n" << endl;
			}
			break;
        	}
	} // flanks map should now be filled with every flank found.
	//unsigned initialFlankNumber = flanks.size();

	if (opt::printFlanks > 0) {
		map<FastaRecord, map<FastaRecord, map<string, int> > >:: iterator read1_it;
		map<FastaRecord, map<string, int> >::iterator read2_it;

		string read1OutputPath(opt::outputPrefix);
		read1OutputPath.append("_reads_1.fq");
		ofstream read1Stream(read1OutputPath.c_str());
		assert_good(read1Stream, read1OutputPath);

		string read2OutputPath(opt::outputPrefix);
		read2OutputPath.append("_reads_2.fq");
		ofstream read2Stream(read2OutputPath.c_str());
		assert_good(read2Stream, read2OutputPath);

		for (read1_it = flanks.begin(); read1_it != flanks.end(); read1_it++) {
			FastaRecord read1 = read1_it->first;
			for (read2_it = flanks[read1].begin(); read2_it != flanks[read1].end(); read2_it++) {
				FastaRecord read2 = read2_it->first;

				read1Stream << ">" << read1.id.substr(0,read2.id.length()-2)
					<< "_" << read2_it->second["startposition"]
					<< "_" << read2_it->second["endposition"]
						- read2_it->second["startposition"] << "/1\n";
				read1Stream << read1.seq << endl;

				read2Stream << ">" << read2.id.substr(0,read2.id.length()-2)
					<< "_" << read2_it->second["startposition"]
					<< "_" << read2_it->second["endposition"]
						- read2_it->second["startposition"] << "/2\n";
				read2Stream << read2.seq << endl;
			}
		}
		assert_good(read1Stream, read1OutputPath.c_str());
		read1Stream.close();
		assert_good(read2Stream, read2OutputPath.c_str());
		read2Stream.close();
	}

//	exit(EXIT_SUCCESS);

	/** map for merged sequence resutls */
	map<string, map<int, map<string, string> > > allmerged;
	unsigned gapsclosed=0;

	for (unsigned i = 0; i<opt::kvector.size(); i++) {
		opt::k = opt::kvector.at(i);
		Kmer::setLength(opt::k);

	//	Counters g_count;
		BloomFilter bloom;

		if (!opt::bloomFilterPaths.empty() && i <= opt::bloomFilterPaths.size()) {

			if (opt::verbose)
				std::cerr << "Loading bloom filter from `"
					<< opt::bloomFilterPaths.at(i) << "'...\n";

			const char* inputPath = opt::bloomFilterPaths.at(i).c_str();
			ifstream inputBloom(inputPath, ios_base::in | ios_base::binary);
			assert_good(inputBloom, inputPath);
			inputBloom >> bloom;
			assert_good(inputBloom, inputPath);
			inputBloom.close();
		} else {
			if (opt::verbose > 0)
				cerr << "Building bloom filter" << endl;
			size_t bits = opt::bloomSize * 8 / 2;
			bloom = BloomFilter(bits);
#ifdef _OPENMP
			ConcurrentBloomFilter<BloomFilter> cbf(bloom, 1000);
			for (int i = optind; i < argc; i++)
				Bloom::loadFile(cbf, opt::k, argv[i], 0 /*opt::verbose*/);
#else
			for (int i = optind; i < argc; i++)
				Bloom::loadFile(bloom, opt::k, argv[i], 0/* opt::verbose*/);
#endif
		}

		DBGBloom<BloomFilter> g(bloom);

		if (opt::verbose > 0)
			cerr << "Starting K run with k = " << opt::k << endl;
		kRun(params, opt::k, g, allmerged, flanks, gapsclosed);

		if (opt::cascade == 1) {
			if (opt::verbose > 0)
				cerr << "Building cascading bloom filter" << endl;
			// Specify bloom filter size in bits.
			size_t bits = opt::bloomSize * 8 / 2;
			CascadingBloomFilter tempBloom(bits);
#ifdef _OPENMP
			ConcurrentBloomFilter<CascadingBloomFilter> cbf(tempBloom, 1000);
			for (int i = optind; i < argc; i++)
				Bloom::loadFile(cbf, opt::k, string(argv[i]), 0/*opt::verbose*/);
#else
			for (int i = optind; i < argc; i++)
				Bloom::loadFile(tempBloom, opt::k, string(argv[i]), 0/*opt::verbose*/);
#endif
			bloom = tempBloom.getBloomFilter(tempBloom.MAX_COUNT-1);

			DBGBloom<BloomFilter> g(bloom);

			if (opt::verbose > 0)
				cerr << "Starting K run with cascading bloom filter" << endl;
			kRun(params, opt::k, g, allmerged, flanks, gapsclosed);
		}
		if (opt::verbose > 0)
			cerr << "k" << opt::k << " run complete\n\n";

	}

	if (opt::verbose > 0) {
		cerr << "K sweep complete" << endl;
		cerr << "Starting new scaffold creation..." << endl;
	}

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
	if (opt::verbose > 0) {
		cerr << "New scaffold complete" << endl;
		cerr << "Gaps closed = " << gapsclosed << endl;
		cerr << (float)100 * gapsclosed / gapsfound << "%\n\n";
	}

//	if (opt::verbose > 0) {
//		cerr <<
//		//	"Processed " << initialFlankNumber << " pseudoreads\n"
//			"Merged (Unique path + Multiple paths): "
//				<< g_count.uniquePath + g_count.multiplePaths
//				<< /*" (" << setprecision(3) <<  (float)100
//				    * (g_count.uniquePath + g_count.multiplePaths) /
//				   initialFlankNumber
//				<< "%)*/"\n"
//			"No start/goal kmer: " << g_count.noStartOrGoalKmer
//				<< /*" (" << setprecision(3) << (float)100
//					* g_count.noStartOrGoalKmer / initialFlankNumber
//				<< "%)*/"\n"
//			"No path: " << g_count.noPath
//				<< /*" (" << setprecision(3) << (float)100
//					* g_count.noPath / initialFlankNumber
//				<< "%)*/"\n"
//			"Unique path: " << g_count.uniquePath
//				<</* " (" << setprecision(3) << (float)100
//					* g_count.uniquePath / initialFlankNumber
//				<< "%)*/"\n"
//			"Multiple paths: " << g_count.multiplePaths
//				<</* " (" << setprecision(3) << (float)100
//					* g_count.multiplePaths / initialFlankNumber
//				<< "%)*/"\n"
//			"Too many paths: " << g_count.tooManyPaths
//				<</* " (" << setprecision(3) << (float)100
//					* g_count.tooManyPaths / initialFlankNumber
//				<< "%)*/"\n"
//			"Too many branches: " << g_count.tooManyBranches
//				<</* " (" << setprecision(3) << (float)100
//					* g_count.tooManyBranches / initialFlankNumber
//				<< "%)*/"\n"
//			"Too many path/path mismatches: " << g_count.tooManyMismatches
//				<</* " (" << setprecision(3) << (float)100
//					* g_count.tooManyMismatches / initialFlankNumber
//				<< "%)*/"\n"
//			"Too many path/read mismatches: " << g_count.tooManyReadMismatches
//				<< /*" (" << setprecision(3) << (float)100
//					* g_count.tooManyReadMismatches / initialFlankNumber
//				<< "%)*/"\n"
//			"Contains cycle: " << g_count.containsCycle
//				<< /*" (" << setprecision(3) << (float)100
//					* g_count.containsCycle / initialFlankNumber
//				<< "%)*/"\n"
//			"Exceeded mem limit: " << g_count.exceededMemLimit
//				<</* " (" << setprecision(3) << (float)100
//					* g_count.exceededMemLimit / initialFlankNumber
//				<< "%)*/"\n"
//			"Skipped: " << g_count.skipped
//				<</* " (" << setprecision(3) << (float)100
//					* g_count.skipped / initialFlankNumber
//				<< "%)*/"\n"
//		/*	"Bloom filter FPR: " << setprecision(3) << 100 * bloom->FPR()
//				<< "%\n"*/;
//	}
	assert_good(scaffoldStream, scaffoldOutputPath.c_str());
	scaffoldStream.close();
	assert_good(mergedStream, mergedOutputPath.c_str());
	mergedStream.close();

	if (!opt::dotPath.empty()) {
		assert_good(dotStream, opt::dotPath);
		dotStream.close();
	}

	if (!opt::tracefilePath.empty()) {
		assert_good(traceStream, opt::tracefilePath);
		traceStream.close();
	}

//	delete bloom;

	return 0;
}
