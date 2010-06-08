#include "Estimate.h"
#include "Histogram.h"
#include "SAM.h"
#include "Stats.h"
#include "Uncompress.h"
#include <cassert>
#include <cerrno>
#include <climits>
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

	/** Reverse-forward mate pair orientation. */
	static bool rf = false;

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
	distanceList.reserve(pairs.size());
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

static void writeEstimate(ostream& out,
		const ContigNode& id0, const ContigNode& id1,
		unsigned len0, unsigned len1,
		const AlignPairVec& pairs, const PDF& pdf)
{
	if (pairs.size() < opt::npairs)
		return;

	Estimate est;
	est.contig = id1;
	est.distance = estimateDistance(len0, len1,
			pairs, pdf, est.numPairs);
	est.stdDev = pdf.getSampleStdDev(est.numPairs);

	if (est.numPairs >= opt::npairs) {
		if (opt::dot) {
			if (id0.sense())
				est.contig.flip();
			out << '"' << id0 << "\" -> " << est << '\n';
		} else
			out << ' ' << est;
	} else if (opt::verbose > 1) {
		cerr << "warning: " << id0 << ',' << id1 << ' '
			<< est.numPairs << " of " << pairs.size()
			<< " pairs fit the expected distribution\n";
	}
}

/** Generate distance estimates for the specified alignments. */
static void writeEstimates(ostream& out,
		const vector<SAMRecord>& pairs,
		const vector<unsigned>& lengthVec, const PDF& pdf)
{
	assert(!pairs.empty());
	LinearNumKey id0 = g_contigIDs.serial(pairs.front().rname);
	assert(id0 < lengthVec.size());
	unsigned len0 = lengthVec[id0];
	if (len0 < opt::seedLen)
		return; // Skip contigs shorter than the seed length.

	if (!opt::dot)
		out << pairs.front().rname;

	typedef map<ContigNode, AlignPairVec> Pairs;
	Pairs dataMap[2];
	for (AlignPairVec::const_iterator it = pairs.begin();
			it != pairs.end(); ++it)
		dataMap[it->isReverse()][ContigNode(it->mrnm,
				it->isReverse() == it->isMateReverse())]
			.push_back(*it);

	for (int sense0 = false; sense0 <= true; sense0++) {
		if (!opt::dot && sense0 == 1)
			out << " ;";
		for (Pairs::const_iterator it = dataMap[sense0].begin();
				it != dataMap[sense0].end(); ++it)
			writeEstimate(out, ContigNode(id0, sense0), it->first,
					len0, lengthVec[it->first.id()],
					it->second, pdf);
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
	unsigned numRF = distanceHist.count(INT_MIN, 0);
	unsigned numFR = distanceHist.count(1, INT_MAX);
	unsigned numTotal = distanceHist.size();
	cout << "Mate orientation FR: " << numFR << setprecision(3)
		<< " (" << (float)100*numFR/numTotal << "%)"
		<< " RF: " << numRF << setprecision(3)
		<< " (" << (float)100*numRF/numTotal << "%)"
		<< endl;
	if (numFR < numRF) {
		cout << "The mate pairs of this library are oriented "
			"reverse-forward (RF)." << endl;
		opt::rf = true;
	}
	assert(!opt::rf);

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
