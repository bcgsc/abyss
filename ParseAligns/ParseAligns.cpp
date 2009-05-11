#include "Aligner.h"
#include "PairUtils.h"
#include "Stats.h"
#include <algorithm>
#include <cerrno>
#include <climits> // for INT_MIN
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <functional>
#include <getopt.h>
#include <iostream>
#include <sstream>
#include <string>

using namespace std;

#define PROGRAM "ParseAligns"

static const char *VERSION_MESSAGE =
PROGRAM " (ABySS) " VERSION "\n"
"Written by Jared Simpson and Shaun Jackman.\n"
"\n"
"Copyright 2009 Canada's Michael Smith Genome Science Centre\n";

static const char *USAGE_MESSAGE =
"Usage: " PROGRAM " [OPTION]... [FILE]...\n"
"Write read pairs that align to the same contig to FRAGMENTS or HISTOGRAM.\n"
"Write read pairs that align to different contigs to standard output.\n"
"Alignments may be in FILE(s) or standard input.\n"
"\n"
"  -k, --kmer=KMER_SIZE  k-mer size\n"
"  -f, --frag=FRAGMENTS  write fragment sizes to this file\n"
"  -h, --hist=HISTOGRAM  write the fragment size histogram to this file\n"
"  -v, --verbose         display verbose output\n"
"      --help            display this help and exit\n"
"      --version         output version information and exit\n"
"\n"
"Report bugs to <" PACKAGE_BUGREPORT ">.\n";

namespace opt {
	static int k;
	static int verbose;
	static string fragPath;
	static string histPath;
}

static const char* shortopts = "k:f:h:v";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
	{ "kmer",    required_argument, NULL, 'k' },
	{ "frag",    required_argument, NULL, 'f' },
	{ "hist",    required_argument, NULL, 'h' },
	{ "verbose", no_argument,       NULL, 'v' },
	{ "help",    no_argument,       NULL, OPT_HELP },
	{ "version", no_argument,       NULL, OPT_VERSION },
	{ NULL, 0, NULL, 0 }
};

static struct {
	int alignments;
	int numDifferent;
	int numSame;
	int numInvalid;
	int numMissed;
	int numMulti;
	int numNonSingle;
} stats;
static ostream& pairedAlignFile = cout;
static ofstream fragFile;
static Histogram histogram;

// TYPEDEFS
typedef hash_map<string, AlignmentVector> ReadAlignMap;

// FUNCTIONS
bool checkUniqueAlignments(int kmer, const AlignmentVector& alignVec);

std::string makePairID(std::string refID);

/**
 * Return the size of the fragment demarcated by the specified
 * alignments.
 * @return INT_MIN if the pair is invalid
 */
static int fragmentSize(const Alignment& a0, const Alignment& a1)
{
	assert(a0.contig == a1.contig);
	if (a0.isRC == a1.isRC)
		return INT_MIN;
	const Alignment& f = a0.isRC ? a1 : a0;
	const Alignment& r = a0.isRC ? a0 : a1;
	return r - f;
}

static void handleAlignmentPair(ReadAlignMap::const_iterator iter,
		ReadAlignMap::const_iterator pairIter,
		const string& currID, const string& pairID)
{
	// Both reads must align to a unique location.
	// The reads are allowed to span more than one contig, but
	// at least one of the two reads must span no more than
	// two contigs.
	bool isRefUnique = checkUniqueAlignments(opt::k,
			iter->second);
	bool isPairUnique = checkUniqueAlignments(opt::k,
			pairIter->second);
	const unsigned MAX_SPAN = 2;
	if ((iter->second.size() <= MAX_SPAN
				|| pairIter->second.size() <= MAX_SPAN)
			&& isRefUnique && isPairUnique) {
		// Iterate over the vectors, outputting the aligments
		for(AlignmentVector::const_iterator refAlignIter = iter->second.begin(); refAlignIter != iter->second.end(); ++refAlignIter)
		{
			for(AlignmentVector::const_iterator pairAlignIter = pairIter->second.begin(); pairAlignIter != pairIter->second.end(); ++pairAlignIter)
			{
				// Are they on the same contig and the ONLY alignments?
				if(refAlignIter->contig == pairAlignIter->contig)
				{
					if((iter->second.size() == 1 && pairIter->second.size() == 1))
					{
						int size = fragmentSize(
								*refAlignIter,
								*pairAlignIter);
						if (size > INT_MIN) {
							histogram.addDataPoint(size);
							if (!opt::fragPath.empty()) {
								fragFile << size << "\n";
								assert(fragFile.good());
							}
							stats.numSame++;
						} else
							stats.numInvalid++;
					}
				}
				else
				{
					// Print the alignment and the swapped alignment
					pairedAlignFile << currID << " " << *refAlignIter << " " << pairID << " " << *pairAlignIter << "\n";
					pairedAlignFile << pairID << " " << *pairAlignIter << " " << currID << " " << *refAlignIter << "\n";
					assert(pairedAlignFile.good());
					stats.numDifferent++;
				}
			}
		}
	}
	else
	{
		if(!isRefUnique || !isPairUnique)
		{
			stats.numMulti++;
		}
		else
		{
			stats.numNonSingle++;
		}
	}
}

// Read in the alignments file into the table
static void readAlignments(istream& in, ReadAlignMap* pout)
{
	ReadAlignMap& out = *pout;
	string line;
	for (string line; getline(in, line);) {
		stringstream s(line);
		string readID;
		s >> readID;
		AlignmentVector& alignments = out[readID];
		if (!alignments.empty()) {
			cerr << "error: duplicate read ID `"
				<< readID << "'\n";
			exit(EXIT_FAILURE);
		}
		for (Alignment ali; s >> ali;)
			alignments.push_back(ali);
		stats.alignments++;

		string pairID = makePairID(readID);

		// Find the pair align
		ReadAlignMap::iterator iter = out.find(readID);
		ReadAlignMap::iterator pairIter = out.find(pairID);
		if(pairIter != out.end()) {
			handleAlignmentPair(iter, pairIter, readID, pairID);

			// Erase the pair as its not needed (faster to mark it as invalid?)
			out.erase(pairIter);
			out.erase(iter);
		}
	}
	assert(in.eof());
}

static void assert_open(ifstream& f, const string& p)
{
	if (!f.is_open()) {
		cerr << p << ": " << strerror(errno) << endl;
		exit(EXIT_FAILURE);
	}
}

static void readAlignmentsFile(string path, ReadAlignMap* pout)
{
	if (opt::verbose > 0)
		cerr << "Reading `" << path << "'..." << endl;
	ifstream fin(path.c_str());
	assert_open(fin, path);
	readAlignments(fin, pout);
	fin.close();
}

int main(int argc, char* const* argv)
{
	bool die = false;
	for (char c; (c = getopt_long(argc, argv,
					shortopts, longopts, NULL)) != -1;) {
		istringstream arg(optarg != NULL ? optarg : "");
		switch (c) {
			case '?': die = true; break;
			case 'k': arg >> opt::k; break;
			case 'f': arg >> opt::fragPath; break;
			case 'h': arg >> opt::histPath; break;
			case 'v': opt::verbose++; break;
			case OPT_HELP:
				cout << USAGE_MESSAGE;
				exit(EXIT_SUCCESS);
			case OPT_VERSION:
				cout << VERSION_MESSAGE;
				exit(EXIT_SUCCESS);
		}
	}

	if (opt::k <= 0) {
		cerr << PROGRAM ": " << "missing -k,--kmer option\n";
		die = true;
	}

	if (opt::fragPath.empty() && opt::histPath.empty()) {
		cerr << PROGRAM ": " << "missing either -f,--frag or "
			"-h,--hist option\n";
		die = true;
	}

	if (die) {
		cerr << "Try `" << PROGRAM
			<< " --help' for more information.\n";
		exit(EXIT_FAILURE);
	}

	if (!opt::fragPath.empty()) {
		fragFile.open(opt::fragPath.c_str());
		assert(fragFile.is_open());
	}

	ReadAlignMap alignTable;
	if (optind < argc) {
		for_each(argv + optind, argv + argc,
				bind2nd(ptr_fun(readAlignmentsFile), &alignTable));
	} else {
		if (opt::verbose > 0)
			cerr << "Reading from standard input..." << endl;
		readAlignments(cin, &alignTable);
	}
	if (opt::verbose > 0)
		cerr << "Read " << stats.alignments << " alignments" << endl;

	stats.numMissed = alignTable.size();

	if (opt::verbose > 0)
		cerr << "Unmatched: " << stats.numMissed
			<< " Same: " << stats.numSame
			<< " Invalid: " << stats.numInvalid
			<< " Diff: " << stats.numDifferent
			<< " Multi: " << stats.numMulti
			<< " Non-singular: " << stats.numNonSingle
			<< endl;

	if (!opt::fragPath.empty())
		fragFile.close();

	if (!opt::histPath.empty()) {
		ofstream histFile(opt::histPath.c_str());
		assert(histFile.is_open());
		histFile << histogram;
		assert(histFile.good());
		histFile.close();
	}
}

bool checkUniqueAlignments(int kmer, const AlignmentVector& alignVec)
{
	// Ensure that no read start position hit to more than 1 contig
	assert(!alignVec.empty());
	
	const int num_starts = alignVec.front().read_length - kmer + 1;
	int* coverage = new int[num_starts];
	
	for(int i = 0; i < num_starts; ++i)
	{
		coverage[i] = 0;
	}
	
	for(AlignmentVector::const_iterator iter = alignVec.begin(); iter != alignVec.end(); ++iter)
	{
		int length = iter->align_length;
		int start = iter->read_start_pos;
		for(int i = 0; i < (length - kmer + 1); ++i)
		{
			int start_idx = start + i;
			assert(start_idx >= 0 && start_idx < num_starts);
			coverage[start_idx]++;
		}
	}
	
	bool unique = true;
	//printf("Coverage: \n");
	for(int i = 0; i < num_starts; ++i)
	{
		//printf("%d", coverage[i]);
		if(coverage[i] > 1)
		{
			unique = false;
		}
	}
	//printf("\n");
	delete [] coverage;
	return unique;
}

// Change the id into the id of its pair
std::string makePairID(std::string refID)
{
	std::string pairID = refID;
	// Change the last character
	ssize_t lastIdx = pairID.size() - 1;
	assert(lastIdx > 0);
	char c = refID[lastIdx];
	bool matched = true;
	switch (c) {
		case '1': c = '2'; break;
		case '2': c = '1'; break;
		case 'A': c = 'B'; break;
		case 'B': c = 'A'; break;
		case 'F': c = 'R'; break;
		case 'R': c = 'F'; break;
		default: matched = false; break;
	}
	if (matched) {
		pairID[lastIdx] = c;
	} else {
		static const char *suffix0 = "forward";
		static const char *suffix1 = "reverse";
		ssize_t i = refID.length() - strlen(suffix0);
		assert(i > 0);
		string suffix = refID.substr(i);
		if (suffix == suffix0) {
			pairID.replace(i, string::npos, suffix1);
		} else if (suffix == suffix1) {
			pairID.replace(i, string::npos, suffix0);
		} else {
			cerr << "error: read ID `" << refID
				<< "' must end in one of\n"
				"\t/1 and /2 or"
				" _A and _B or"
				" _F and _R or"
				" _forward and _reverse\n";
			exit(EXIT_FAILURE);
		}
	}
	return pairID;
}
