#include "Estimate.h"
#include "Histogram.h"
#include "PairUtils.h"
#include "SAM.h"
#include "Stats.h"
#include "Uncompress.h"
#include <cassert>
#include <cerrno>
#include <cstdlib>
#include <cstring> // for strerror
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <limits> // for numeric_limits
#include <sstream>
#include <string>
#include <vector>

using namespace std;

#define PROGRAM "DistanceEst"

static const char VERSION_MESSAGE[] =
PROGRAM " (" PACKAGE_NAME ") " VERSION "\n"
"Written by Jared Simpson and Shaun Jackman.\n"
"\n"
"Copyright 2010 Canada's Michael Smith Genome Science Centre\n";

static const char USAGE_MESSAGE[] =
"Usage: " PROGRAM " [OPTION]... HIST [PAIR]\n"
"Estimate distances between contigs using paired-end alignments.\n"
"  HIST  distribution of fragments size\n"
"  PAIR  alignments between contigs\n"
"\n"
"  -k, --kmer=KMER_SIZE  k-mer size\n"
"  -n, --npairs=NPAIRS   minimum number of pairs\n"
"  -s, --seed-length=L   minimum length of the seed contigs [100]\n"
"  -o, --out=FILE        write result to FILE\n"
"      --dot             output overlaps in dot format\n"
"  -v, --verbose         display verbose output\n"
"      --help            display this help and exit\n"
"      --version         output version information and exit\n"
"\n"
"Report bugs to <" PACKAGE_BUGREPORT ">.\n";

namespace opt {
	static unsigned k;

	/** Output in dot format. */
	int dot; // used by Estimate

	static unsigned npairs;
	static unsigned seedLen = 100;
	static int verbose;
	static string out;
}

static const char shortopts[] = "k:n:o:s:v";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
	{ "dot",         no_argument,       &opt::dot, 1, },
	{ "kmer",        required_argument, NULL, 'k' },
	{ "npairs",      required_argument, NULL, 'n' },
	{ "out",         required_argument, NULL, 'o' },
	{ "seed-length", required_argument, NULL, 's' },
	{ "verbose",     no_argument,       NULL, 'v' },
	{ "help",        no_argument,       NULL, OPT_HELP },
	{ "version",     no_argument,       NULL, OPT_VERSION },
	{ NULL, 0, NULL, 0 }
};

typedef vector<SAMRecord> AlignPairVec;

struct PairedData
{
	AlignPairVec pairVec[2];
};

typedef map<ContigID, PairedData> PairDataMap;

/** Estimate the distance between two contigs.
 * @param numPairs [out] the number of pairs that agree with the
 * expected distribution
 * @return the estimated distance
 */
static int estimateDistance(unsigned len0, unsigned len1,
		const AlignPairVec& pairs, const PDF& pdf,
		unsigned& numPairs)
{
	// The provisional fragment sizes are calculated as if the contigs
	// were perfectly adjacent with no overlap or gap.
	vector<int> distanceList;
	for (AlignPairVec::const_iterator it = pairs.begin();
			it != pairs.end(); ++it) {
		Alignment a0 = *it;
		if (a0.isRC)
			a0 = a0.flipTarget(len0);
		int a1 = it->isMateReverse() ? it->mpos : len1 - it->mpos;
		int distance = len0 + a1 - a0.targetAtQueryStart();
		assert(distance > 0);
		distanceList.push_back(distance);
	}
	return maxLikelihoodEst(-opt::k+1, pdf.getMaxIdx(),
			distanceList, pdf, numPairs);
}

static void writeEstimates(ostream& out,
		const vector<SAMRecord>& currPairs,
		const vector<unsigned>& lengthVec, const PDF& pdf)
{
	assert(!currPairs.empty());
	ContigID refContigID = currPairs.front().rname;

	// From this point all ids will be interpreted as integers
	// They must be strictly > 0 and contiguous
	LinearNumKey refNumericID
		= convertContigIDToLinearNumKey(refContigID);
	assert(refNumericID < lengthVec.size());

	// Only process contigs that are a reasonable length
	unsigned refLength = lengthVec[refNumericID];
	if (refLength < opt::seedLen)
		return;

	if (!opt::dot)
		out << refContigID;

	// Seperate the pairings by direction (pairs aligning in the
	// same comp as the contig are sense pairs) and by the contig
	// they align to
	for(size_t dirIdx = 0; dirIdx <= 1; ++dirIdx)
	{
		// If this is the second direction, write a seperator
		if (!opt::dot && dirIdx == 1)
			out << " ;";
		ContigNode refContig(refContigID, dirIdx);

		PairDataMap dataMap;
		for (AlignPairVec::const_iterator iter = currPairs.begin();
				iter != currPairs.end(); ++iter) {
			if (iter->isReverse() == (bool)dirIdx) {
				PairedData& pd = dataMap[iter->mrnm];
				size_t compIdx = (size_t)iter->isMateReverse();
				assert(compIdx < 2);
				pd.pairVec[compIdx].push_back(*iter);
			}
		}

		// For each contig that is paired, compute the distance
		for (PairDataMap::const_iterator pdIter = dataMap.begin();
				pdIter != dataMap.end(); ++pdIter) {
			const ContigID& pairID = pdIter->first;
			// Check if the pairs are in a valid orientation
			if (pdIter->second.pairVec[0].size() >= opt::npairs
					&& pdIter->second.pairVec[1].size()
					>= opt::npairs) {
				cerr << "warning: inconsistent pairing between "
					<< refContig << ' '
					<< pairID << '+' << ' '
					<< pdIter->second.pairVec[1].size()
					<< ' '
					<< pairID << '-' << ' '
					<< pdIter->second.pairVec[0].size()
					<< '\n';
				continue;
			}

			unsigned pairDirIdx = pdIter->second.pairVec[0].size()
				>= opt::npairs ? 0 : 1;
			const AlignPairVec& pairVec
				= pdIter->second.pairVec[pairDirIdx];
			unsigned numPairs = pairVec.size();
			if (numPairs >= opt::npairs) {
				Estimate est;
				est.contig = ContigNode(
						pairID, dirIdx == pairDirIdx);
				est.distance = estimateDistance(
						refLength, lengthVec.at(est.contig.id()),
						pairVec, pdf, est.numPairs);
				est.stdDev = pdf.getSampleStdDev(est.numPairs);

				if (est.numPairs >= opt::npairs) {
					if (opt::dot) {
						if (dirIdx)
							est.contig.flip();
						out << '"' << refContig << "\" -> "
							<< est << '\n';
					} else
						out << ' ' << est;
				} else if (opt::verbose > 1) {
					cerr << "warning: "
						<< refContigID << (dirIdx ? '-' : '+')
						<< ','
						<< est.contig << ' '
						<< est.numPairs << " of "
						<< numPairs << " pairs"
						" fit the expected distribution\n";
				}
			}
		}
	}
	if (!opt::dot)
		out << "\n";
	assert(out.good());
}

static void assert_open(ifstream& f, const string& p)
{
	if (f.is_open())
		return;
	cerr << p << ": " << strerror(errno) << endl;
	exit(EXIT_FAILURE);
}

/** Load a histogram from the specified file. */
static Histogram loadHist(const string& path)
{
	ifstream in(path.c_str());
	assert_open(in, path);

	Histogram hist;
	in >> hist;
	assert(in.eof());

	if (hist.empty()) {
		cerr << "error: the histogram `" << path << "' is empty\n";
		exit(EXIT_FAILURE);
	}
	return hist;
}

/** Read contig lengths from SAM headers. */
static void readContigLengths(istream& in, vector<unsigned>& lengths)
{
	assert(in);
	assert(lengths.empty());
	assert(g_contigIDs.empty());
	for (string line; in.peek() == '@' && getline(in, line);) {
		istringstream ss(line);
		string type, tag;
		ss >> type;
		if (type != "@SQ")
			continue;
		ss >> ws;

		getline(ss, tag, ':');
		assert(tag == "SN");
		string id;
		ss >> id >> ws;

		getline(ss, tag, ':');
		assert(tag == "LN");
		unsigned len;
		ss >> len;

		assert(ss);
		(void)g_contigIDs.serial(id);
		lengths.push_back(len);
	}
	assert(!lengths.empty());
}

int main(int argc, char** argv)
{
	bool die = false;
	for (int c; (c = getopt_long(argc, argv,
					shortopts, longopts, NULL)) != -1;) {
		istringstream arg(optarg != NULL ? optarg : "");
		switch (c) {
			case '?': die = true; break;
			case 'k': arg >> opt::k; break;
			case 'n': arg >> opt::npairs; break;
			case 'o': arg >> opt::out; break;
			case 's': arg >> opt::seedLen; break;
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
		cerr << PROGRAM ": missing -k,--kmer option\n";
		die = true;
	}

	if (opt::npairs <= 0) {
		cerr << PROGRAM ": missing -n,--npairs option\n";
		die = true;
	}

	if (argc - optind < 1) {
		cerr << PROGRAM ": missing arguments\n";
		die = true;
	} else if (argc - optind > 2) {
		cerr << PROGRAM ": too many arguments\n";
		die = true;
	}

	if (die) {
		cerr << "Try `" << PROGRAM
			<< " --help' for more information.\n";
		exit(EXIT_FAILURE);
	}

	if (opt::seedLen < 2*opt::k)
		cerr << "warning: the seed-length should be at least twice k:"
			" k=" << opt::k << ", s=" << opt::seedLen << '\n';

	string distanceCountFile(argv[optind++]);
	string alignFile(argv[optind] == NULL ? "-" : argv[optind++]);

	ifstream inFile(alignFile.c_str());
	istream& in(strcmp(alignFile.c_str(), "-") == 0 ? cin : inFile);

	if (strcmp(alignFile.c_str(), "-") != 0)
		assert_open(inFile, alignFile);

	ofstream outFile;
	if (!opt::out.empty()) {
		outFile.open(opt::out.c_str());
		assert(outFile.is_open());
	}
	ostream& out = opt::out.empty() ? cout : outFile;

	if (opt::dot)
		out << "digraph dist {\n"
			"k=" << opt::k << "\t"
			"n=" << opt::npairs << "\t"
			"s=" << opt::seedLen << "\n";

	// Read the contig lengths.
	vector<unsigned> contigLens;
	readContigLengths(in, contigLens);
	g_contigIDs.lock();

	// Read the fragment size distribution.
	Histogram distanceHist = loadHist(distanceCountFile);
	distanceHist.eraseNegative();
	Histogram trimmedHist = distanceHist.trimFraction(0.0001);
	PDF empiricalPDF(trimmedHist);

	// Estimate the distances between contigs.
	vector<SAMRecord> alignments(1);
	in >> alignments.front();
	assert(in);
	for (SAMRecord sam; in >> sam; alignments.push_back(sam)) {
		if (sam.rname != alignments.front().rname) {
			writeEstimates(out, alignments, contigLens, empiricalPDF);
			alignments.clear();
		}
	}
	assert(in.eof());
	writeEstimates(out, alignments, contigLens, empiricalPDF);

	if (opt::dot)
		out << "}\n";

	inFile.close();
	if (!opt::out.empty())
		outFile.close();

	return 0;
}
