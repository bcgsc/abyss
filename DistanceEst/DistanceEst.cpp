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
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <getopt.h>
#include <iomanip>
#include <iostream>
#include <limits> // for numeric_limits
#include <vector>
#if _OPENMP
# include <omp.h>
#endif
#include "DataBase/Options.h"
#include "DataBase/DB.h"

using namespace std;

#define PROGRAM "DistanceEst"

DB db;

static const char VERSION_MESSAGE[] =
PROGRAM " (" PACKAGE_NAME ") " VERSION "\n"
"Written by Jared Simpson and Shaun Jackman.\n"
"\n"
"Copyright 2014 Canada's Michael Smith Genome Sciences Centre\n";

static const char USAGE_MESSAGE[] =
"Usage: " PROGRAM " -k<kmer> -s<seed-length> -n<npairs> [OPTION]... HIST [PAIR]\n"
"Estimate distances between contigs using paired-end alignments.\n"
"\n"
" Arguments:\n"
"\n"
"  HIST  distribution of fragments size\n"
"  PAIR  alignments between contigs\n"
"\n"
" Options:\n"
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
"      --mle             use the MLE [default]\n"
"                        (maximum likelihood estimator)\n"
"      --median          use the difference of the population median\n"
"                        and the sample median\n"
"      --mean            use the difference of the population mean\n"
"                        and the sample mean\n"
"      --dist            output the graph in dist format [default]\n"
"      --dot             output the graph in GraphViz format\n"
"      --gv              output the graph in GraphViz format\n"
"      --gfa             output the graph in GFA2 format\n"
"      --gfa2            output the graph in GFA2 format\n"
"  -j, --threads=N       use N parallel threads [1]\n"
"  -v, --verbose         display verbose output\n"
"      --help            display this help and exit\n"
"      --version         output version information and exit\n"
"      --db=FILE         specify path of database repository in FILE\n"
"      --library=NAME    specify library NAME for sqlite\n"
"      --strain=NAME     specify strain NAME for sqlite\n"
"      --species=NAME    specify species NAME for sqlite\n"
"\n"
"Report bugs to <" PACKAGE_BUGREPORT ">.\n";

/** Which estimator to use. See opt::method. */
enum { MLE, MEAN, MEDIAN };

namespace opt {
	string db;
	dbVars metaVars;
	unsigned k; // used by Estimate.h

	/** Output graph format. */
	int format = DIST;

	/** Minimum distance between contigs. */
	static int minDist = numeric_limits<int>::min();

	/** Maximum distance between contigs. */
	static int maxDist = numeric_limits<int>::max();

	static unsigned seedLen;
	static unsigned npairs;
	static unsigned minMapQ = 10;

	/** Reverse-forward mate pair orientation. */
	static int rf = -1;

	/** Which estimator to use. */
	static int method = MLE;

	static int verbose;
	static string out;
	static int threads = 1;
}

static const char shortopts[] = "j:k:l:n:o:q:s:v";

enum { OPT_HELP = 1, OPT_VERSION,
	OPT_MIND, OPT_MAXD, OPT_FR, OPT_RF,
	OPT_DB, OPT_LIBRARY, OPT_STRAIN, OPT_SPECIES
};
//enum { OPT_HELP = 1, OPT_VERSION,
//	OPT_MIND, OPT_MAXD, OPT_FR, OPT_RF
//};

static const struct option longopts[] = {
	{ "dist",        no_argument,       &opt::format, DIST, },
	{ "dot",         no_argument,       &opt::format, DOT, },
	{ "gv",          no_argument,       &opt::format, DOT, },
	{ "gfa",         no_argument,       &opt::format, GFA2, },
	{ "gfa2",        no_argument,       &opt::format, GFA2, },
	{ "fr",          no_argument,       &opt::rf, false },
	{ "rf",          no_argument,       &opt::rf, true },
	{ "min-align",   required_argument, NULL, 'l' },
	{ "mind",        required_argument, NULL, OPT_MIND },
	{ "maxd",        required_argument, NULL, OPT_MAXD },
	{ "mle",         no_argument,       &opt::method, MLE },
	{ "median",      no_argument,       &opt::method, MEDIAN },
	{ "mean",        no_argument,       &opt::method, MEAN },
	{ "kmer",        required_argument, NULL, 'k' },
	{ "npairs",      required_argument, NULL, 'n' },
	{ "out",         required_argument, NULL, 'o' },
	{ "min-mapq",    required_argument, NULL, 'q' },
	{ "seed-length", required_argument, NULL, 's' },
	{ "threads",     required_argument, NULL, 'j' },
	{ "verbose",     no_argument,       NULL, 'v' },
	{ "help",        no_argument,       NULL, OPT_HELP },
	{ "version",     no_argument,       NULL, OPT_VERSION },
	{ "db",          required_argument, NULL, OPT_DB },
	{ "library",     required_argument, NULL, OPT_LIBRARY },
	{ "strain",      required_argument, NULL, OPT_STRAIN },
	{ "species",     required_argument, NULL, OPT_SPECIES },
	{ NULL, 0, NULL, 0 }
};

/** A collection of aligned read pairs. */
typedef vector<SAMRecord> Pairs;

/** Estimate the distance between two contigs using the difference of
 * the population mean and the sample mean.
 * @param numPairs [out] the number of pairs that agree with the
 * expected distribution
 * @return the estimated distance
 */
static int estimateDistanceUsingMean(
		const std::vector<int>& samples, const PMF& pmf,
		unsigned& numPairs)
{
	Histogram h(samples.begin(), samples.end());
	int d = (int)round(pmf.mean() - h.mean());

	// Count the number of samples that agree with the distribution.
	unsigned n = 0;
	for (Histogram::const_iterator it = h.begin();
			it != h.end(); ++it)
		if (pmf[it->first + d] > pmf.minProbability())
			n += it->second;

	numPairs = n;
	return d;
}

/** Estimate the distance between two contigs using the difference of
 * the population median and the sample median.
 * @param numPairs [out] the number of pairs that agree with the
 * expected distribution
 * @return the estimated distance
 */
static int estimateDistanceUsingMedian(
		const std::vector<int>& samples, const PMF& pmf,
		unsigned& numPairs)
{
	Histogram h(samples.begin(), samples.end());
	int d = (int)round(pmf.median() - h.median());
	// Count the number of samples that agree with the distribution.
	unsigned n = 0;
	for (Histogram::const_iterator it = h.begin();
			it != h.end(); ++it)
		if (pmf[it->first + d] > pmf.minProbability())
			n += it->second;

	numPairs = n;
	return d;
}

/** Global variable to track a recommended minAlign parameter */
unsigned g_recMA;

static struct {
	/* Fragment stats are considered only for fragments aligning
	 * to different contigs, and where the contig is >=opt::seedLen. */
	unsigned total_frags;
	unsigned dup_frags;
} stats;

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
	unsigned orig = fragments.size();
	sort(fragments.begin(), fragments.end());
	fragments.erase(unique(fragments.begin(), fragments.end()),
			fragments.end());
	numPairs = fragments.size();
	assert((int)orig - (int)numPairs >= 0);
	stats.total_frags += orig;
	stats.dup_frags += orig - numPairs;

	if (numPairs < opt::npairs)
		return INT_MIN;

	vector<int> fragmentSizes;
	fragmentSizes.reserve(fragments.size());
	unsigned ma = opt::minAlign;
	for (Fragments::const_iterator it = fragments.begin();
			it != fragments.end(); ++it) {
		int x = it->second - it->first;
		if (!opt::rf && opt::method == MLE
				&& x <= 2 * int(ma - 1)) {
			unsigned align = x / 2;
			if (opt::verbose > 0)
#pragma omp critical(cerr)
				cerr << PROGRAM ": warning: The observed fragment of "
					"size " << x << " bp is shorter than 2*l "
					"(l=" << opt::minAlign << ").\n";
			ma = min(ma, align);
		}
		fragmentSizes.push_back(x);
	}

#pragma omp critical(g_recMA)
	g_recMA = min(g_recMA, ma);
	switch (opt::method) {
	  case MLE:
		// Use the maximum likelihood estimator.
		return maximumLikelihoodEstimate(ma,
				opt::minDist, opt::maxDist,
				fragmentSizes, pmf, len0, len1, opt::rf, numPairs);
	  case MEAN:
		// Use the difference of the population mean
		// and the sample mean.
		return estimateDistanceUsingMean(
				fragmentSizes, pmf, numPairs);
	  case MEDIAN:
		// Use the difference of the population median
		// and the sample median.
		return estimateDistanceUsingMedian(
				fragmentSizes, pmf, numPairs);
	  default:
		assert(false);
		abort();
	}
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
		} else if (opt::format == GFA2) {
			// Output only one of the two complementary edges.
			if (len1 < opt::seedLen || e.first < e.second || e.first == e.second)
#pragma omp critical(out)
				out << "G\t*"
					<< '\t' << get(g_contigNames, e.first)
					<< '\t' << get(g_contigNames, e.second)
					<< '\t' << est.distance
					<< '\t' << (int)ceilf(est.stdDev)
					<< "\tFC:i:" << est.numPairs
					<< '\n';
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
	ContigID id0(get(g_contigNames, pairs.front().rname));
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
		dataMap[it->isReverse()][find_vertex(
				it->mrnm, it->isReverse() == it->isMateReverse(),
				g_contigNames)]
			.push_back(*it);

	for (int sense0 = false; sense0 <= true; sense0++) {
		if (opt::format == DIST && sense0)
			ss << " ;";
		const PairsMap& x = dataMap[sense0 ^ opt::rf];
		for (PairsMap::const_iterator it = x.begin();
				it != x.end(); ++it)
			writeEstimate(opt::format == DIST ? ss : out,
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
			&& get(g_contigNames, it->rname)
				< get(g_contigNames, out.front().rname)) {
		cerr << "error: input must be sorted: saw `"
			<< out.front().rname << "' before `"
			<< it->rname << "'\n";
		exit(EXIT_FAILURE);
	}
}

int main(int argc, char** argv)
{
	if (!opt::db.empty())
		opt::metaVars.resize(3);

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
			case OPT_DB:
				arg >> opt::db; break;
			case OPT_LIBRARY:
				arg >> opt::metaVars[0]; break;
			case OPT_STRAIN:
				arg >> opt::metaVars[1]; break;
			case OPT_SPECIES:
				arg >> opt::metaVars[2]; break;
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

	if (!opt::db.empty()) {
		init(db,
				opt::db,
				opt::verbose,
				PROGRAM,
				opt::getCommand(argc, argv),
				opt::metaVars
		);
	}

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
	else if (opt::format == GFA2)
		out << "H\tVN:Z:2.0\n";

	vector<int> vals = make_vector<int>()
		<< opt::k
		<< opt::seedLen
		<< opt::npairs;

	vector<string> keys = make_vector<string>()
		<< "K"
		<< "SeedLen"
		<< "NumPairs";

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

	vals += make_vector<int>()
		<< numFR
		<< numRF;

	keys += make_vector<string>()
		<< "FR_orientation"
		<< "RF_orientation";

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

	vals += make_vector<int>()
		<< opt::minDist
		<< opt::maxDist
		<< (int)round(h.mean())
		<< h.median()
		<< (int)round(h.sd())
		<< h.size()
		<< h.minimum()
		<< h.maximum();

	keys += make_vector<string>()
		<< "minDist"
		<< "maxDist"
		<< "mean"
		<< "median"
		<< "sd"
		<< "n"
		<< "min"
		<< "max";

	// Read the contig lengths.
	vector<unsigned> contigLens;

	vals += make_vector<int>()
		<< readContigLengths(in, contigLens);

	keys += make_vector<string>()
		<< "CntgCounted";

//	readContigLengths(in, contigLens);

	g_contigNames.lock();

	// Estimate the distances between contigs.
	istream_iterator<SAMRecord> it(in), last;
	if (contigLens.size() == 1) {
		// When mapping to a single contig, no alignments spanning
		// contigs are expected.
		assert(in.eof());
		exit(EXIT_SUCCESS);
	}
	assert(in);

	g_recMA = opt::minAlign;
#pragma omp parallel
	for (vector<SAMRecord> records;;) {
		records.clear();
#pragma omp critical(in)
		readPairs(it, last, records);
		if (records.empty())
			break;
		writeEstimates(out, records, contigLens, pmf);
	}

	if (opt::verbose > 0) {
		float prop_dups = (float)100 * stats.dup_frags / stats.total_frags;
		cerr << "Duplicate rate of spanning fragments: "
			<< stats.dup_frags << "/"
			<< stats.total_frags << " ("
			<< setprecision(3) << prop_dups << "%)\n";
		if (prop_dups > 50)
			cerr << PROGRAM << ": warning: duplicate rate of fragments "
				"spanning more than one contig is high.\n";
	}

	vals += make_vector<int>()
		<< stats.total_frags
		<< stats.dup_frags;

		keys += make_vector<string>()
			<< "total_frags"
			<< "dupl_frags";

	if (!opt::db.empty()) {
		for (unsigned i=0; i<vals.size(); i++)
			addToDb(db, keys[i], vals[i]);
	}

	if (opt::verbose > 0 && g_recMA != opt::minAlign)
		cerr << PROGRAM << ": warning: MLE will be more accurate if "
			"l is decreased to " << g_recMA << ".\n";

	assert(in.eof());

	if (opt::format == DOT)
		out << "}\n";
	return 0;
}
