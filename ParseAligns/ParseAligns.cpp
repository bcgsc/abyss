#include "Aligner.h"
#include "Histogram.h"
#include "PairUtils.h"
#include "SAM.h"
#include "Uncompress.h"
#include <algorithm>
#include <cerrno>
#include <climits> // for INT_MIN
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <functional>
#include <getopt.h>
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

#define PROGRAM "ParseAligns"

static const char VERSION_MESSAGE[] =
PROGRAM " (" PACKAGE_NAME ") " VERSION "\n"
"Written by Jared Simpson and Shaun Jackman.\n"
"\n"
"Copyright 2009 Canada's Michael Smith Genome Science Centre\n";

static const char USAGE_MESSAGE[] =
"Usage: " PROGRAM " [OPTION]... [FILE]...\n"
"Write read pairs that align to the same contig to FRAGMENTS or HISTOGRAM.\n"
"Write read pairs that align to different contigs to standard output.\n"
"Alignments may be in FILE(s) or standard input.\n"
"\n"
"  -k, --kmer=KMER_SIZE  k-mer size\n"
"  -d, --dist=DISTANCE   write distance estimates to this file\n"
"  -f, --frag=FRAGMENTS  write fragment sizes to this file\n"
"  -h, --hist=HISTOGRAM  write the fragment size histogram to this file\n"
"      --sam             alignments are in SAM format\n"
"      --kaligner        alignments are in KAligner format\n"
"  -c, --cover=COVERAGE  coverage cut-off for distance estimates\n"
"  -v, --verbose         display verbose output\n"
"      --help            display this help and exit\n"
"      --version         output version information and exit\n"
"\n"
"Report bugs to <" PACKAGE_BUGREPORT ">.\n";

namespace opt {
	static int k;
	static unsigned c;
	static int verbose;
	static string distPath;
	static string fragPath;
	static string histPath;

	/** Input alignment format. */
	static int inputFormat;
	enum { KALIGNER, SAM };
}

static const char shortopts[] = "d:k:f:h:c:v";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
	{ "dist",    required_argument, NULL, 'd' },
	{ "kmer",    required_argument, NULL, 'k' },
	{ "frag",    required_argument, NULL, 'f' },
	{ "hist",    required_argument, NULL, 'h' },
	{ "kaligner",no_argument, &opt::inputFormat, opt::KALIGNER },
	{ "sam",     no_argument, &opt::inputFormat, opt::SAM },
	{ "cover",   required_argument, NULL, 'c' },
	{ "verbose", no_argument,       NULL, 'v' },
	{ "help",    no_argument,       NULL, OPT_HELP },
	{ "version", no_argument,       NULL, OPT_VERSION },
	{ NULL, 0, NULL, 0 }
};

static struct {
	int alignments;
	int numDifferent;
	int numSame;
	int numMisoriented;
	int numMissed;
	int numMulti;
	int numSplit;
} stats;

static ostream& pairedAlignFile = cout;
static ofstream fragFile;
static Histogram histogram;

// TYPEDEFS
typedef hash_map<string, AlignmentVector> ReadAlignMap;
typedef hash_map<ContigID, EstimateRecord> EstimateMap;

static EstimateMap estMap;

static bool checkUniqueAlignments(const AlignmentVector& alignVec);
string makePairID(string id);

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

static void addEstimate(EstimateMap& map, const Alignment& a, Estimate& est, bool reverse)
{
	//count up the number of estimates that agree
	bool placed = false;
	bool a_isRC = a.isRC != reverse;
	EstimateMap::iterator estimatesIt = map.find(a.contig);
	if (estimatesIt != map.end()) {
		EstimateVector& estimates =
			estimatesIt->second.estimates[a_isRC];
		for (EstimateVector::iterator estIt = estimates.begin();
				estIt != estimates.end(); ++estIt) {
			if (estIt->nID == est.nID) {
				estIt->numPairs++;
				estIt->distance += est.distance;
				placed = true;
				break;
			}
		}
	}
	if (!placed)
		map[a.contig].estimates[a_isRC].push_back(est);

}

static void doReadIntegrity(const ReadAlignMap::value_type& a)
{
	AlignmentVector::const_iterator refAlignIter = a.second.begin();
	unsigned firstStart, lastEnd, largestSize;
	Alignment first, last, largest;

	firstStart = refAlignIter->read_start_pos;
	lastEnd = firstStart + refAlignIter->align_length;
	largestSize = refAlignIter->align_length;
	first = last = largest = *refAlignIter;
	++refAlignIter;

	//for each alignment in the vector a.second
	for (; refAlignIter != a.second.end(); ++refAlignIter) {
		if ((unsigned)refAlignIter->read_start_pos < firstStart) {
			firstStart = refAlignIter->read_start_pos;
			first = *refAlignIter;
		}
		if ((unsigned)(refAlignIter->read_start_pos +
					refAlignIter->align_length) > lastEnd) {
			lastEnd = refAlignIter->read_start_pos +
				refAlignIter->align_length;
			last = *refAlignIter;
		}
		if ((unsigned)refAlignIter->align_length > largestSize) {
			largestSize = refAlignIter->align_length;
			largest = *refAlignIter;
		}
	}

	if (largest.contig != last.contig) {
		Estimate est;
		unsigned largest_end = 
			largest.read_start_pos + largest.align_length - opt::k;
		int distance = last.read_start_pos - largest_end;
		est.nID = convertContigIDToLinearNumKey(last.contig);
		est.distance = distance - opt::k;
		est.numPairs = 1;
		est.stdDev = 0;
		//weird file format...
		est.isRC = largest.isRC != last.isRC;

		addEstimate(estMap, largest, est, false);
	}

	if (largest.contig != first.contig) {
		largest.flipQuery();
		first.flipQuery();
		Estimate est;
		unsigned largest_end = 
			largest.read_start_pos + largest.align_length - opt::k;
		int distance = first.read_start_pos - largest_end;
		est.nID = convertContigIDToLinearNumKey(first.contig);
		est.distance = distance - opt::k;
		est.numPairs = 1;
		est.stdDev = 0;
		//weird file format...
		est.isRC = largest.isRC != first.isRC;

		addEstimate(estMap, largest, est, false);
	}


#if 0
	//for each alignment in the vector a.second
	for (AlignmentVector::const_iterator refAlignIter = a.second.begin();
			refAlignIter != a.second.end(); ++refAlignIter) {
		//for each alignment after the current one
		for (AlignmentVector::const_iterator alignIter = a.second.begin();
				alignIter != a.second.end(); ++alignIter) {
			//make sure both alignments aren't for the same contig
			if (alignIter->contig != refAlignIter->contig) {
				Estimate est;
				//Make sure the distance is read as 0 if the two contigs are
				//directly adjacent to each other. A -ve number suggests an
				//overlap.
				assert(refAlignIter->read_start_pos != alignIter->read_start_pos);
				const Alignment& a = refAlignIter->read_start_pos < alignIter->read_start_pos ? *refAlignIter : *alignIter;
				const Alignment& b = refAlignIter->read_start_pos > alignIter->read_start_pos ? *refAlignIter : *alignIter;
				unsigned a_end = a.read_start_pos + a.align_length - opt::k;
				int distance = b.read_start_pos - a_end;
				est.nID = convertContigIDToLinearNumKey(b.contig);
				est.distance = distance - opt::k;
				est.numPairs = 1;
				est.stdDev = 0;
				//weird file format...
				est.isRC = a.isRC != b.isRC;

				addEstimate(estMap, a, est, false);
			}
		}
	}
#endif
}

static void generateDistFile()
{
	ofstream distFile(opt::distPath.c_str());
	assert(distFile.is_open());
	for (EstimateMap::iterator mapIt = estMap.begin();
			mapIt != estMap.end(); ++mapIt) {
		//Skip empty iterators
		assert(!mapIt->second.estimates[0].empty() || !mapIt->second.estimates[1].empty());
		distFile << mapIt->first << " :";
		for (int refIsRC = 0; refIsRC <= 1; refIsRC++) {
			if (refIsRC)
				distFile << " |";

			for (EstimateVector::iterator vecIt = mapIt->second.estimates[refIsRC].begin();
					vecIt != mapIt->second.estimates[refIsRC].end(); ++vecIt) {
				vecIt->distance = (int)round((double)vecIt->distance /
						(double)vecIt->numPairs);
				if (vecIt->numPairs >= opt::c && vecIt->numPairs != 0
						/*&& vecIt->distance > 1 - opt::k*/)
					distFile << ' ' << *vecIt;
			}
		}
		distFile << '\n';
		assert(distFile.good());
	}
	distFile.close();
}

static bool needsFlipping(const string& id);

/**
 * Return an alignment flipped as necessary to produce an alignment
 * pair whose expected orientation is forward-reverse.  If the
 * expected orientation is forward-forward, then reverse the first
 * alignment, so that the alignment is forward-reverse, which is
 * required by DistanceEst.
 */
static const Alignment flipAlignment(const Alignment& a,
		const string& id)
{
	return needsFlipping(id) ? a.flipQuery() : a;
}

static void handleAlignmentPair(const ReadAlignMap::value_type& curr,
		const ReadAlignMap::value_type& pair)
{
	const string& currID = curr.first;
	const string& pairID = pair.first;

	// Both reads must align to a unique location.
	// The reads are allowed to span more than one contig, but
	// at least one of the two reads must span no more than
	// two contigs.
	const unsigned MAX_SPAN = 2;
	if (curr.second.empty() || pair.second.empty()) {
		stats.numMissed++;
	} else if (!checkUniqueAlignments(curr.second)
			|| !checkUniqueAlignments(pair.second)) {
		stats.numMulti++;
	} else if (curr.second.size() > MAX_SPAN
			&& pair.second.size() > MAX_SPAN) {
		stats.numSplit++;
	} else {
		// Iterate over the vectors, outputting the aligments
		bool counted = false;
		for (AlignmentVector::const_iterator refAlignIter
					= curr.second.begin();
				refAlignIter != curr.second.end(); ++refAlignIter) {
			for (AlignmentVector::const_iterator pairAlignIter
						= pair.second.begin();
					pairAlignIter != pair.second.end();
					++pairAlignIter) {
				const Alignment& a0 = flipAlignment(*refAlignIter,
						currID);
				const Alignment& a1 = flipAlignment(*pairAlignIter,
						pairID);

				// Are they on the same contig and the ONLY alignments?
				if (a0.contig == a1.contig) {
					int size = fragmentSize(a0, a1);
					if (curr.second.size() == 1
							&& pair.second.size() == 1) {
						if (size > INT_MIN) {
							histogram.insert(size);
							if (!opt::fragPath.empty()) {
								fragFile << size << "\n";
								assert(fragFile.good());
							}
							stats.numSame++;
						} else
							stats.numMisoriented++;
						counted = true;
					}
					if (size < opt::k) {
						/* For an inverted repeat, print only the
						 * one alignment because only one distance
						 * estimate is expected. For a non-inverted
						 * repeat, two distance esimates are expected
						 * (both to and from the repeat to itself) so
						 * print both alignments. */
						pairedAlignFile
							<< currID << ' ' << a0 << ' '
							<< pairID << ' ' << a1 << '\n';
						if (a0.isRC != a1.isRC)
							pairedAlignFile
								<< pairID << ' ' << a1 << ' '
								<< currID << ' ' << a0 << '\n';
						assert(pairedAlignFile.good());
					}
				} else {
					// Print the alignment and the swapped alignment
					pairedAlignFile
						<< currID << ' ' << a0 << ' '
						<< pairID << ' ' << a1 << '\n'
						<< pairID << ' ' << a1 << ' '
						<< currID << ' ' << a0 << '\n';
					assert(pairedAlignFile.good());
				}
			}
		}
		if (!counted)
			stats.numDifferent++;
	}
}

static void printStatus(const ReadAlignMap& map)
{
	if (opt::verbose == 0)
		return;

	static size_t prevBuckets;
	size_t buckets = map.bucket_count();
	if (stats.alignments % 1000000 == 0 || buckets != prevBuckets) {
		prevBuckets = buckets;
		size_t size = map.size();
		cerr << "Read " << stats.alignments << " alignments. "
			<< "Hash load: " << size << " / " << buckets
			<< " = " << (float)size / buckets << endl;
	}
}

static void handleAlignment(
		const ReadAlignMap::value_type& alignments,
		ReadAlignMap& out)
{
	stats.alignments++;

	string pairID = makePairID(alignments.first);
	ReadAlignMap::iterator pairIter = out.find(pairID);
	if (pairIter != out.end()) {
		handleAlignmentPair(*pairIter, alignments);
		out.erase(pairIter);
	} else if (!out.insert(alignments).second) {
		cerr << "error: duplicate read ID `" << alignments.first
			<< "'\n";
		exit(EXIT_FAILURE);
	}

	if (!opt::distPath.empty() && alignments.second.size() >= 2)
		doReadIntegrity(alignments);

	printStatus(out);
}

// Read in the alignments file into the table
static void readAlignments(istream& in, ReadAlignMap* pout)
{
	string line;
	for (string line; getline(in, line);) {
		istringstream s(line);
		switch (opt::inputFormat) {
		  case opt::SAM:
		  {
			SAMRecord sam;
			s >> sam;
			assert(s);
			ReadAlignMap::value_type alignments(
					sam.qname, AlignmentVector());
			if (!(sam.flag & SAMRecord::FUNMAP))
				alignments.second.push_back(sam);
			handleAlignment(alignments, *pout);
			break;
		  }
		  case opt::KALIGNER:
		  {
			string readID;
			s >> readID;
			ReadAlignMap::value_type alignments(
					readID, AlignmentVector());
			copy(istream_iterator<Alignment>(s),
					istream_iterator<Alignment>(),
					back_inserter(alignments.second));
			handleAlignment(alignments, *pout);
			break;
		  }
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
	for (int c; (c = getopt_long(argc, argv,
					shortopts, longopts, NULL)) != -1;) {
		istringstream arg(optarg != NULL ? optarg : "");
		switch (c) {
			case '?': die = true; break;
			case 'k': arg >> opt::k; break;
			case 'c': arg >> opt::c; break;
			case 'd': arg >> opt::distPath; break;
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

	if (opt::verbose > 0)
		cerr << "Mateless: " << alignTable.size()
			<< " Unaligned: " << stats.numMissed
			<< " Same: " << stats.numSame
			<< " Misoriented: " << stats.numMisoriented
			<< " Diff: " << stats.numDifferent
			<< " Multi: " << stats.numMulti
			<< " Split: " << stats.numSplit
			<< endl;

	if (!opt::distPath.empty())
		generateDistFile();

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

static bool checkUniqueAlignments(const AlignmentVector& alignVec)
{
	// Ensure that no read start position hit to more than 1 contig
	assert(!alignVec.empty());

	const int num_starts = alignVec.front().read_length - opt::k + 1;
	vector<unsigned> coverage(num_starts);

	for(AlignmentVector::const_iterator iter = alignVec.begin(); iter != alignVec.end(); ++iter)
	{
		int length = iter->align_length;
		int start = iter->read_start_pos;
		for (int i = 0; i < length - opt::k + 1; i++) {
			int start_idx = start + i;
			assert(start_idx >= 0 && start_idx < num_starts);
			coverage[start_idx]++;
		}
	}

	bool unique = true;
	for(int i = 0; i < num_starts; ++i)
	{
		if(coverage[i] > 1)
			unique = false;
	}
	return unique;
}

static bool endsWith(const string& s, const string& suffix)
{
	ssize_t i = s.length() - suffix.length();
	return i < 0 ? false : s.substr(i) == suffix;
}

static bool replaceSuffix(string& s,
		const string& suffix0, const string& suffix1)
{
	if (endsWith(s, suffix0)) {
		s.replace(s.length() - suffix0.length(), string::npos,
				suffix1);
		return true;
	} else if (endsWith(s, suffix1)) {
		s.replace(s.length() - suffix1.length(), string::npos,
				suffix0);
		return true;
	} else
		return false;
}

/** Return the mate ID of the specified read ID. */
string makePairID(string id)
{
	assert(!id.empty());
	char& c = id[id.length() - 1];
	switch (c) {
		case '1': c = '2'; return id;
		case '2': c = '1'; return id;
		case 'A': c = 'B'; return id;
		case 'B': c = 'A'; return id;
		case 'F': c = 'R'; return id;
		case 'R': c = 'F'; return id;
		case 'f': c = 'r'; return id;
		case 'r': c = 'f'; return id;
	}

	if (replaceSuffix(id, "forward", "reverse")
				|| replaceSuffix(id, "F3", "R3"))
		return id;

	cerr << "error: read ID `" << id << "' must end in one of\n"
		"\t1 and 2 or A and B or F and R or"
		" F3 and R3 or forward and reverse\n";
	exit(EXIT_FAILURE);
}

static bool needsFlipping(const string& id)
{
	return endsWith(id, "F3");
}
