#include "Estimate.h"
#include "Histogram.h"
#include "IOUtil.h"
#include "MLE.h"
#include "PMF.h"
#include "SAM.h"
#include "Uncompress.h"
#include "Graph/Options.h" // for opt::k
#include <algorithm>
#include <cassert>
#include <climits>
#include <cstdlib>
#include <fstream>
#include <getopt.h>
#include <iomanip>
#include <iostream>
#include <iterator> // for istream_iterator
#include <limits> // for numeric_limits
#include <sstream>
#include <string>
#include <vector>
#if _OPENMP
# include <omp.h>
#endif

using namespace std;

#define PROGRAM "DistanceEst"

static const char VERSION_MESSAGE[] =
PROGRAM " (" PACKAGE_NAME ") " VERSION "\n"
"Written by Jared Simpson and Shaun Jackman.\n"
"\n"
"Copyright 2012 Canada's Michael Smith Genome Science Centre\n";

static const char USAGE_MESSAGE[] =
"Usage: " PROGRAM " [OPTION]... HIST [PAIR]\n"
"Estimate distances between contigs using paired-end alignments.\n"
"  HIST  distribution of fragments size\n"
"  PAIR  alignments between contigs\n"
"\n"
"      --mind=N          minimum distance between contigs [-(k-1)]\n"
"      --maxd=N          maximum distance between contigs\n"
"      --fr              force the orientation to forward-reverse\n"
"      --rf              force the orientation to reverse-forward\n"
"  -k, --kmer=N          set --mind to -(k-1) bp\n"
"  -l, --min-align=N     the minimal alignment size [1]\n"
"  -n, --npairs=NPAIRS   minimum number of pairs\n"
"  -s, --seed-length=L   minimum length of the seed contigs\n"
"  -q, --min-mapq=N      ignore alignments with mapping quality\n"
"                        less than this threshold [10]\n"
"  -o, --out=FILE        write result to FILE\n"
"      --dist            output graph in dist format [default]\n"
"      --dot             output graph in dot format\n"
"  -j, --threads=N       use N parallel threads [1]\n"
"  -v, --verbose         display verbose output\n"
"      --help            display this help and exit\n"
"      --version         output version information and exit\n"
"\n"
"Report bugs to <" PACKAGE_BUGREPORT ">.\n";

namespace opt {
	unsigned k; // used by Estimate.h

	/** Output graph format. */
	int format = DIST;

	/** The minimal alignment size. */
	static int minAlign = 1;

	/** Minimum distance between contigs. */
	static int minDist = numeric_limits<int>::min();

	/** Maximum distance between contigs. */
	static int maxDist = numeric_limits<int>::max();

	static unsigned seedLen;
	static unsigned npairs;
	static unsigned minMapQ = 10;

	/** Reverse-forward mate pair orientation. */
	static int rf = -1;

	static int verbose;
	static string out;
	static int threads = 1;
}

static const char shortopts[] = "j:k:l:n:o:q:s:v";

enum { OPT_HELP = 1, OPT_VERSION,
	OPT_MIND, OPT_MAXD, OPT_FR, OPT_RF
};

static const struct option longopts[] = {
	{ "dist",        no_argument,       &opt::format, DIST, },
	{ "dot",         no_argument,       &opt::format, DOT, },
	{ "fr",          no_argument,       &opt::rf, false },
	{ "rf",          no_argument,       &opt::rf, true },
	{ "min-align",   required_argument, NULL, 'l' },
	{ "mind",        required_argument, NULL, OPT_MIND },
	{ "maxd",        required_argument, NULL, OPT_MAXD },
	{ "kmer",        required_argument, NULL, 'k' },
	{ "npairs",      required_argument, NULL, 'n' },
	{ "out",         required_argument, NULL, 'o' },
	{ "min-mapq",    required_argument, NULL, 'q' },
	{ "seed-length", required_argument, NULL, 's' },
	{ "threads",     required_argument,	NULL, 'j' },
	{ "verbose",     no_argument,       NULL, 'v' },
	{ "help",        no_argument,       NULL, OPT_HELP },
	{ "version",     no_argument,       NULL, OPT_VERSION },
	{ NULL, 0, NULL, 0 }
};

/** A collection of aligned read pairs. */
typedef vector<SAMRecord> Pairs;

/** Estimate the distance between two contigs.
 * @param numPairs [out] the number of pairs that agree with the
 * expected distribution
 * @return the estimated distance
 */
static int estimateDistance(unsigned len0, unsigned len1,
		const Pairs& pairs, const PMF& pmf,
		unsigned& numPairs)
{
	// The provisional fragment sizes are calculated as if the contigs
	// were perfectly adjacent with no overlap or gap.
	typedef vector<pair<int, int> > Fragments;
	Fragments fragments;
	fragments.reserve(pairs.size());
	for (Pairs::const_iterator it = pairs.begin();
			it != pairs.end(); ++it) {
		int a0 = it->targetAtQueryStart();
		int a1 = it->mateTargetAtQueryStart();
		if (it->isReverse())
			a0 = len0 - a0;
		if (!it->isMateReverse())
			a1 = len1 - a1;
		fragments.push_back(opt::rf
				? make_pair(a1, len1 + a0)
				: make_pair(a0, len0 + a1));
	}

	// Remove duplicate fragments.
	sort(fragments.begin(), fragments.end());
	fragments.erase(unique(fragments.begin(), fragments.end()),
			fragments.end());
	numPairs = fragments.size();
	if (numPairs < opt::npairs)
		return INT_MIN;

	vector<int> fragmentSizes;
	fragmentSizes.reserve(fragments.size());
	for (Fragments::const_iterator it = fragments.begin();
			it != fragments.end(); ++it) {
		int x = it->second - it->first;
		if (!opt::rf && x <= 2 * int(opt::minAlign - 1)) {
			cerr << PROGRAM ": error: The observed fragment of size "
				<< x << " bp is shorter than 2*l "
				"(l=" << opt::minAlign << "). "
				"Decrease l to " << x / 2 << ".\n";
			exit(EXIT_FAILURE);
		}
		fragmentSizes.push_back(x);
	}

	return maximumLikelihoodEstimate(opt::minAlign,
			opt::minDist, opt::maxDist,
			fragmentSizes, pmf, len0, len1, opt::rf, numPairs);
}

static void writeEstimate(ostream& out,
		const ContigNode& id0, const ContigNode& id1,
		unsigned len0, unsigned len1,
		const Pairs& pairs, const PMF& pmf)
{
	if (pairs.size() < opt::npairs)
		return;

	DistanceEst est;
	est.distance = estimateDistance(len0, len1,
			pairs, pmf, est.numPairs);
	est.stdDev = pmf.getSampleStdDev(est.numPairs);

	std::pair<ContigNode, ContigNode> e(id0, id1 ^ id0.sense());
	if (est.numPairs >= opt::npairs) {
		if (opt::format == DOT) {
#pragma omp critical(out)
			out << get(g_contigNames, e) << " [" << est << "]\n";
		} else
			out << ' ' << get(g_contigNames, id1) << ',' << est;
	} else if (opt::verbose > 1) {
#pragma omp critical(cerr)
		cerr << "warning: " << get(g_contigNames, e)
			<< " [d=" << est.distance << "] "
			<< est.numPairs << " of " << pairs.size()
			<< " pairs fit the expected distribution\n";
	}
}

/** Generate distance estimates for the specified alignments. */
static void writeEstimates(ostream& out,
		const vector<SAMRecord>& pairs,
		const vector<unsigned>& lengthVec, const PMF& pmf)
{
	assert(!pairs.empty());
	ContigID id0(pairs.front().rname);
	assert(id0 < lengthVec.size());
	unsigned len0 = lengthVec[id0];
	if (len0 < opt::seedLen)
		return; // Skip contigs shorter than the seed length.

	ostringstream ss;
	if (opt::format == DIST)
		ss << pairs.front().rname;

	typedef map<ContigNode, Pairs> PairsMap;
	PairsMap dataMap[2];
	for (Pairs::const_iterator it = pairs.begin();
			it != pairs.end(); ++it)
		dataMap[it->isReverse()][ContigNode(it->mrnm,
				it->isReverse() == it->isMateReverse())]
			.push_back(*it);

	for (int sense0 = false; sense0 <= true; sense0++) {
		if (opt::format == DIST && sense0)
			ss << " ;";
		const PairsMap& x = dataMap[sense0 ^ opt::rf];
		for (PairsMap::const_iterator it = x.begin();
				it != x.end(); ++it)
			writeEstimate(opt::format == DOT ? out : ss,
					ContigNode(id0, sense0), it->first,
					len0, lengthVec[it->first.id()],
					it->second, pmf);
	}
	if (opt::format == DIST)
#pragma omp critical(out)
		out << ss.str() << '\n';
	assert(out.good());
}

/** Load a histogram from the specified file. */
static Histogram loadHist(const string& path)
{
	ifstream in(path.c_str());
	assert_good(in, path);

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
	assert(g_contigNames.empty());
	for (string line; in.peek() == '@' && getline(in, line);) {
		istringstream ss(line);
		string type;
		ss >> type;
		if (type != "@SQ")
			continue;

		string s;
		unsigned len;
		ss >> expect(" SN:") >> s >> expect(" LN:") >> len;
		assert(ss);

		put(g_contigNames, lengths.size(), s);
		lengths.push_back(len);
	}
	if (lengths.empty()) {
		cerr << PROGRAM ": error: no @SQ records in the SAM header\n";
		exit(EXIT_FAILURE);
	}
}

/** Copy records from [it, last) to out and stop before alignments to
 * the next target sequence.
 * @param[in,out] it an input iterator
 */
template<typename It>
static void readPairs(It& it, const It& last, vector<SAMRecord>& out)
{
	assert(out.empty());
	for (; it != last; ++it) {
		if (it->isUnmapped() || it->isMateUnmapped()
				|| !it->isPaired() || it->rname == it->mrnm
				|| it->mapq < opt::minMapQ)
			continue;
		if (!out.empty() && out.back().rname != it->rname)
			break;

		out.push_back(*it);
		SAMRecord& sam = out.back();
		// Clear unused fields.
		sam.qname.clear();
#if SAM_SEQ_QUAL
		sam.seq.clear();
		sam.qual.clear();
#endif
	}

	// Check that the input is sorted.
	if (it != last && !out.empty()
			&& ContigID(it->rname) < ContigID(out.front().rname)) {
		cerr << "error: input must be sorted: saw `"
			<< out.front().rname << "' before `"
			<< it->rname << "'\n";
		exit(EXIT_FAILURE);
	}
}

int main(int argc, char** argv)
{
	bool die = false;
	for (int c; (c = getopt_long(argc, argv,
					shortopts, longopts, NULL)) != -1;) {
		istringstream arg(optarg != NULL ? optarg : "");
		switch (c) {
			case '?': die = true; break;
			case OPT_MIND:
				arg >> opt::minDist;
				break;
			case OPT_MAXD:
				arg >> opt::maxDist;
				break;
			case 'l':
				arg >> opt::minAlign;
				break;
			case 'j': arg >> opt::threads; break;
			case 'k': arg >> opt::k; break;
			case 'n': arg >> opt::npairs; break;
			case 'o': arg >> opt::out; break;
			case 'q': arg >> opt::minMapQ; break;
			case 's': arg >> opt::seedLen; break;
			case 'v': opt::verbose++; break;
			case OPT_HELP:
				cout << USAGE_MESSAGE;
				exit(EXIT_SUCCESS);
			case OPT_VERSION:
				cout << VERSION_MESSAGE;
				exit(EXIT_SUCCESS);
		}
		if (optarg != NULL && !arg.eof()) {
			cerr << PROGRAM ": invalid option: `-"
				<< (char)c << optarg << "'\n";
			exit(EXIT_FAILURE);
		}
	}

	if (opt::k <= 0) {
		cerr << PROGRAM ": missing -k,--kmer option\n";
		die = true;
	}

	if (opt::seedLen <= 0) {
		cerr << PROGRAM ": missing -s,--seed-length option\n";
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

	assert(opt::minAlign > 0);

#if _OPENMP
	if (opt::threads > 0)
		omp_set_num_threads(opt::threads);
#endif

	string distanceCountFile(argv[optind++]);
	string alignFile(argv[optind] == NULL ? "-" : argv[optind++]);

	ifstream inFile(alignFile.c_str());
	istream& in(strcmp(alignFile.c_str(), "-") == 0 ? cin : inFile);

	if (strcmp(alignFile.c_str(), "-") != 0)
		assert_good(inFile, alignFile);

	ofstream outFile;
	if (!opt::out.empty()) {
		outFile.open(opt::out.c_str());
		assert(outFile.is_open());
	}
	ostream& out = opt::out.empty() ? cout : outFile;

	if (opt::format == DOT)
		out << "digraph dist {\ngraph ["
			"k=" << opt::k << " "
			"s=" << opt::seedLen << " "
			"n=" << opt::npairs << "]\n";

	// The fragment size histogram may not be written out until after
	// the alignments complete. Wait for the alignments to complete.
	in.peek();

	// Read the fragment size distribution.
	Histogram distanceHist = loadHist(distanceCountFile);
	unsigned numRF = distanceHist.count(INT_MIN, 0);
	unsigned numFR = distanceHist.count(1, INT_MAX);
	unsigned numTotal = distanceHist.size();
	bool libRF = numFR < numRF;
	if (opt::verbose > 0) {
		cerr << "Mate orientation FR: " << numFR << setprecision(3)
			<< " (" << (float)100*numFR/numTotal << "%)"
			<< " RF: " << numRF << setprecision(3)
			<< " (" << (float)100*numRF/numTotal << "%)\n"
			<< "The library " << distanceCountFile << " is oriented "
			<< (libRF
					? "reverse-forward (RF)" : "forward-reverse (FR)")
			<< ".\n";
	}

	// Determine the orientation of the library.
	if (opt::rf == -1)
		opt::rf = libRF;
	if (opt::rf)
		distanceHist = distanceHist.negate();
	if (opt::rf != libRF)
		cerr << "warning: The orientation is forced to "
			<< (opt::rf
					? "reverse-forward (RF)" : "forward-reverse (FR)")
			<< " which differs from the detected orientation.\n";

	distanceHist.eraseNegative();
	distanceHist.removeNoise();
	distanceHist.removeOutliers();
	Histogram h = distanceHist.trimFraction(0.0001);
	if (opt::verbose > 0)
		cerr << "Stats mean: " << setprecision(4) << h.mean() << " "
			"median: " << setprecision(4) << h.median() << " "
			"sd: " << setprecision(4) << h.sd() << " "
			"n: " << h.size() << " "
			"min: " << h.minimum() << " max: " << h.maximum() << '\n'
			<< h.barplot() << endl;
	PMF pmf(h);

	if (opt::minDist == numeric_limits<int>::min())
		opt::minDist = -opt::k + 1;
	if (opt::maxDist == numeric_limits<int>::max())
		opt::maxDist = pmf.maxValue();
	if (opt::verbose > 0)
		cerr << "Minimum and maximum distance are set to "
			<< opt::minDist << " and " << opt::maxDist << " bp.\n";
	assert(opt::minDist < opt::maxDist);

	// Read the contig lengths.
	vector<unsigned> contigLens;
	readContigLengths(in, contigLens);
	g_contigNames.lock();

	// Estimate the distances between contigs.
	istream_iterator<SAMRecord> it(in), last;
	assert(in);
#pragma omp parallel
	for (vector<SAMRecord> records;;) {
		records.clear();
#pragma omp critical(in)
		readPairs(it, last, records);
		if (records.empty())
			break;
		writeEstimates(out, records, contigLens, pmf);
	}
	assert(in.eof());

	if (opt::format == DOT)
		out << "}\n";
	return 0;
}
