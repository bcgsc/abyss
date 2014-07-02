/**
 * Connect pairs using a Bloom filter de Bruijn graph
 * Copyright 2014 Canada's Michael Smith Genome Science Centre
 */

#include "config.h"

#include "connectpairs/connectpairs.h"
#include "connectpairs/DBGBloom.h"
#include "connectpairs/DBGBloomAlgorithms.h"

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

#define PROGRAM "abyss-connectpairs"

static const char VERSION_MESSAGE[] =
PROGRAM " (" PACKAGE_NAME ") " VERSION "\n"
"Written by Shaun Jackman, Hamid Mohamadi, Anthony Raymond,\n"
"Ben Vandervalk and Daniel Paulino\n"
"\n"
"Copyright 2014 Canada's Michael Smith Genome Science Centre\n";
static const char USAGE_MESSAGE[] =
"Usage: " PROGRAM " -k <kmer size> [-k N ]... -o <output_prefix> -S <path to scaffold file> [options]... <reads1> [reads2]...\n"
"Connect the pairs READS1 and READS2 and close the gap using\n"
"a Bloom filter de Bruijn graph.\n"
"\n"
" Options:\n"
"\n"
"  -S  --input-scaffold=FILE  		load scaffold from FILE\n"
"  -L  --flank-length=N	      		length of flanks to be used as pseudoreads [100]\n"
"  -D  --flank-distance=N     		distance of flank from gap [0]\n"
"  -j, --threads=N            		use N parallel threads [1]\n"
"  -k, --kmer:input-bloom=N:FILE       	the size of a k-mer\n"
"  -b, --bloom-size=N        		size of bloom filter [500M]\n"
"  -B, --max-branches=N    		max branches in de Bruijn graph traversal;\n"
"                             		use 'nolimit' for no limit [350]\n"
"  -d, --dot-file=FILE        		write graph traversals to a DOT file\n"
"  -e, --fix-errors           		find and fix single-base errors when reads\n"
"                             		have no kmers in bloom filter [disabled]\n"
"  -f, --min-frag=N           		min fragment size in base pairs [0]\n"
"  -F, --max-frag=N           		max fragment size in base pairs [1000]\n"
"  -i, --input-bloom=FILE     		load bloom filter from FILE\n"
"  -I, --interleaved          		input reads files are interleaved\n"
"      --mask                 		mask new and changed bases as lower case\n"
"      --no-mask              		do not mask bases [default]\n"
"      --chastity             		discard unchaste reads [default]\n"
"      --no-chastity          		do not discard unchaste reads\n"
"      --trim-masked          		trim masked bases from the ends of reads\n"
"      --no-trim-masked       		do not trim masked bases from the ends\n"
"                             		of reads [default]\n"
"  -l, --long-search          		start path search as close as possible\n"
"                             		to the beginnings of reads. Takes more time\n"
"                             		but improves results when bloom filter false\n"
"                             		positive rate is high [disabled]\n"
"  -m, --read-mismatches=N    		max mismatches between paths and reads; use\n"
"                             		'nolimit' for no limit [nolimit]\n"
"  -M, --max-mismatches=N     		max mismatches between all alternate paths;\n"
"                             		use 'nolimit' for no limit [2]\n"
"  -n  --no-limits            		disable all limits; equivalent to\n"
"                             		'-B nolimit -m nolimit -M nolimit -P nolimit'\n"
"  -o, --output-prefix=FILE   		prefix of output FASTA files [required]\n"
"  -P, --max-paths=N          		merge at most N alternate paths; use 'nolimit'\n"
"                             		for no limit [2]\n"
"  -q, --trim-quality=N       		trim bases from the ends of reads whose\n"
"                             		quality is less than the threshold\n"
"      --standard-quality     		zero quality is `!' (33)\n"
"                             		default for FASTQ and SAM files\n"
"      --illumina-quality     		zero quality is `@' (64)\n"
"                             		default for qseq and export files\n"
"  -r, --read-name=STR        		only process reads with names that contain STR\n"
"  -s, --search-mem=N         		mem limit for graph searches; multiply by the\n"
"                             		number of threads (-j) to get the total mem used\n"
"                             		for graph traversal [500M]\n"
"  -t, --trace-file=FILE      		write graph search stats to FILE\n"
"  -v, --verbose              		display verbose output\n"
"      --help                 		display this help and exit\n"
"      --version              		output version information and exit\n"
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

	/** The size of kmer */
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
} g_count;

static const char shortopts[] = "S:L:D:b:B:d:ef:F:i:Ij:k:lm:M:no:P:q:r:s:t:v";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
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

/** Connect a read pair. */
static void connectPair(const DBGBloom& g,
	const FastqRecord& read1,
	const FastqRecord& read2,
	const ConnectPairsParams& params,
	ofstream& mergedStream,
	ofstream& read1Stream,
	ofstream& read2Stream,
	ofstream& traceStream)
{
	bool skip = false;

#pragma omp atomic
	++g_count.readPairsProcessed;

	if (!opt::readName.empty() &&
		read1.id.find(opt::readName) == string::npos) {
#pragma omp atomic
		++g_count.skipped;
		skip = true;
	}

	if (!skip) {
		/*
		****************************************
		FastqRecord flank1;
		flank1.seq = flank1seq;
		FastqRecord flank2;
		flank2.seq = flank2seq;
		***************************************
		*/
		ConnectPairsResult result = 
			connectPairs(opt::k, read1, read2, g, params);
		
		vector<FastaRecord>& paths = result.mergedSeqs;

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
#pragma omp critical(read1Stream)
					read1Stream << read1;
#pragma omp critical(read2Stream)
					read2Stream << read2;
				}
				else if (paths.size() > 1) {
#pragma omp atomic
					++g_count.multiplePaths;
#pragma omp critical(mergedStream)
					mergedStream << result.consensusSeq;
				}
				else {
#pragma omp atomic
					++g_count.uniquePath;
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
#pragma omp critical(read1Stream)
			read1Stream << read1;
#pragma omp critical(read2Stream)
			read2Stream << read2;
		}
	}

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
				<< "too many path/path mismatches: " << g_count.tooManyMismatches << ", "
				<< "too many path/read mismatches: " << g_count.tooManyReadMismatches << ", "
				<< "contains cycle: " << g_count.containsCycle << ", "
				<< "exceeded mem limit: " << g_count.exceededMemLimit << ", "
				<< "skipped: " << g_count.skipped
				<< ")\n";
		}
	}

}

/** Connect read pairs. */
template <typename FastaStream>
static void connectPairs(const DBGBloom& g,
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
		if (good)
			connectPair(g, a, b, params, mergedStream, read1Stream,
				read2Stream, traceStream);
		else
			break;
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

BloomFilterBase* makeBF (string inputBloomPath) {
	cerr << "makebf begins" << endl;
	BloomFilterBase* bloom = NULL;

	if (!inputBloomPath.empty()) {

		if (opt::verbose)
			std::cerr << "Loading bloom filter from `"
				<< inputBloomPath << "'...\n";

		cout << "before const char" << endl;
		cout << inputBloomPath <<  endl;
		const char* inputPath = inputBloomPath.c_str();
		cout << inputPath << endl;
		cout << "before input bloom" << endl;
		ifstream inputBloom(inputPath, ios_base::in | ios_base::binary);
		assert_good(inputBloom, inputPath);
		cout << "before loaded bloom" << endl;
		BloomFilter* loadedBloom = new BloomFilter();
		cout << "before inputBloom >> loaded bloom" << endl;
		inputBloom >> *loadedBloom;
		cout << "before assert inputBloom, inputPath" << endl;
		assert_good(inputBloom, inputPath);
		cout << "after asser" << endl;
		inputBloom.close();
		cout << "before bloom=loadedBloom" << endl;
		bloom = loadedBloom;
	}
	cout << "makebf exits" <<  endl;
	return bloom;
}

// returns merged sequence resulting from Konnector
string merge(
	const DBGBloom& g,
	unsigned k,
	FastaRecord read1,
	FastaRecord read2,
	const ConnectPairsParams& params)
{
	ConnectPairsResult result = connectPairs(k, read1, read2, g, params);
			
	if (result.pathResult == FOUND_PATH) {
		if (result.mergedSeqs.size() > 1) 
			return result.consensusSeq;
		else
			return result.mergedSeqs.front();	
	}
	else
		return "";
}

void extract(
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
	read1.id = record.id /* + "_" + IntToString(startposition) + "_" + IntToString(endposition - startposition)*/ + "/1";
	// extract right flank
	read2.id = record.id /* + "_" + IntToString(startposition) + "_" + IntToString(endposition - startposition)*/ + "/2";
	string rightflank = seq.substr(endposition, flanklength);
	transform(rightflank.begin(), rightflank.end(), rightflank.begin(), ::toupper);
	read2.seq = reversecompliment(rightflank);
}
// processes one line, finds all gaps of line and runs connector on the gaps. stores results in map. Should pass in a map. This map must be permanently changeable by filler. Return void. 
// TODO: completely redo so that each BF is only loaded once for the entire k sweep
void filler(
	FastaRecord record, 
	int flanklength, 
	const DBGBloom& g, 
	const ConnectPairsParams& params, 
	unsigned k, 
	map<string, map<int, map<string, string> > > &allmerged,
	unsigned &gapnumber) 
{
        int offset = 0;
	int endposition = 0;
       	int startposition = 0;
        string seq = record.seq; //double check if this works	
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
		else if (seq.substr(startposition - flanklength, flanklength).string::find_first_of("Nn") != string::npos
				or seq.substr(endposition, flanklength).string::find_first_of("Nn") != string::npos) {
			offset = endposition;
			continue; 
		}
		
		FastaRecord read1, read2;
		extract(read1, read2, startposition, endposition, seq, flanklength, record);		
		
		cerr << record.id << " found flanks" << endl;
		cerr << read1.seq << endl;
		
		string tempSeq;
		tempSeq = merge(g, k, read1, read2, params);
	
		if (!tempSeq.empty()) {
			allmerged[record.id][startposition]["gap"] = IntToString(endposition - startposition);
			allmerged[record.id][startposition]["seq"] = tempSeq;
			//scaffoldStream << seq.substr(offset, startposition-flanklength-offset) << tempSeq;
			offset = endposition;
		}
		else 	{
			//scaffoldStream << seq.substr(offset, endposition - offset); 
			offset = endposition; //double check
		}
	} 
	//scaffoldStream << seq.substr(offset) << endl;
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
		  case 'i':
			arg >> opt::inputBloomPath; break;
		  case 'I':
			opt::interleaved = true; break;
		  case 'j':
			arg >> opt::threads; break;
		  case 'k': {
			// TODO: double check
			//string tempPath;
			unsigned tempK;
			//char x;
			//cout << "test";
			arg >> tempK; // >> x >> tempPath;
			//cout << ". after arg input" << endl;
			opt::kvector.push_back(tempK);
			//opt::bloomFilterPaths.push_back(tempPath);
			opt::k = tempK;
			break;}
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
		cerr << "creating bloom filter with k value: " << opt::k << endl;
		size_t bits = opt::bloomSize * 8 / 2;
		bloom = new CountingBloomFilter(bits);
		for (int i = optind; i < argc; i++)
			bloom->loadFile(opt::k, string(argv[i]), opt::verbose);

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

	DBGBloom g(*bloom);
	
	string scaffoldOutputPath("newscaffold.fa");	
	ofstream scaffoldStream(scaffoldOutputPath.c_str());
	assert_good(scaffoldStream, scaffoldOutputPath);
	
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
	
	
	// remove fastaInterleave
	/*
	if (opt::interleaved) {
		FastaConcat in(argv + optind, argv + argc,
				FastaReader::FOLD_CASE);
		connectPairs(g, in, params, mergedStream, read1Stream,
				read2Stream, traceStream);
		assert(in.eof());
	} else {
		FastaInterleave in(argv + optind, argv + argc,
				FastaReader::FOLD_CASE);
		connectPairs(g, in, params, mergedStream, read1Stream,
				read2Stream, traceStream);
		assert(in.eof());
	}
	*/

	/** map for merged sequence resutls */
	map<string, map<int, map<string, string> > > allmerged;

	const char* scaffoldInputPath = opt::inputScaffold.c_str();
	FastaReader reader(scaffoldInputPath, FastaReader::FOLD_CASE);
	
	unsigned gapsfound = 0;
	
	// iterates over each read
	cerr << "Entering for loop" << endl;
	for (unsigned i = 0; i<opt::kvector.size(); i++) {
		opt::k = opt::kvector.at(i);
	/*
		
                BloomFilterBase* bloom = NULL;
		if (opt::inputBloomPath.empty()) {	
			cerr << "creating bloom filter with k value: " << opt::k << endl;
			size_t bits = opt::bloomSize * 8 / 2;
			cout << "before bloom = new.." << endl;
                	bloom = new CountingBloomFilter(bits);
			cout << "after bloom = new.." << endl;
                	for (int j = optind; j < argc; j++)
                        	bloom->loadFile(opt::k, string(argv[j]), opt::verbose);
		}
		DBGBloom g(*bloom);	
		 	*/
		cerr << "processing records" << endl;
		for (FastaRecord record;;) {
                	bool good;
                	good = reader >> record;
                	if (good) {
				filler(record, opt::flankLength, g, params, opt::kvector.at(i), allmerged, gapsfound); 
			}
                	else
                        	break;
        	}
		
		
	}
	cerr << "K sweep complete." << endl; 
	cerr << gapsfound << " gaps found." << endl;
	cerr << "Starting new scaffold creation" << endl;
	
	map<string, map<int, map<string, string> > >::iterator scaf_it;
	map<int, map<string, string> >::reverse_iterator pos_it;	
	FastaReader reader2(scaffoldInputPath, FastaReader::FOLD_CASE);
	unsigned gapsclosed = 0;
	
	/** creating new scaffold with gaps closed */
	for (FastaRecord record;;) {
               	bool good;
               	good = reader2 >> record;
               	if (good) {
			scaf_it = allmerged.find(record.id);
			scaffoldStream << record.id << endl;
			string modifiedSeq = record.seq;
			if (scaf_it != allmerged.end()) {
				for (pos_it = allmerged[record.id].rbegin(); pos_it != allmerged[record.id].rend(); pos_it++) {
					modifiedSeq.replace(
						pos_it->first - opt::flankLength, 
						StringToInt(pos_it->second["gap"]) + (opt::flankLength * 2), 
						pos_it->second["seq"]
					);
					cerr << record.id << " gap closed"  << endl;
					cerr << pos_it->second["seq"] << endl;
					gapsclosed++;
				}
				scaffoldStream << modifiedSeq << endl;
			}
			else { 
				scaffoldStream << record.seq  << endl;
			}
		}
               	else
                       	break ;
       	}
	cerr << "New scaffold complete" << endl;
	cerr << gapsclosed << " gaps closed" << endl;
	


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
			"Exceeded mem limit: " << g_count.exceededMemLimit
				<< " (" << setprecision(3) << (float)100
					* g_count.exceededMemLimit / g_count.readPairsProcessed
				<< "%)\n"
			"Skipped: " << g_count.skipped
				<< " (" << setprecision(3) << (float)100
					* g_count.skipped / g_count.readPairsProcessed
				<< "%)\n"
			"Bloom filter FPR: " << setprecision(3) << 100 * bloom->FPR()
				<< "%\n";
	}
	assert_good(scaffoldStream, scaffoldOutputPath.c_str());
	scaffoldStream.close(); 
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

	delete bloom;

	return 0;
}
