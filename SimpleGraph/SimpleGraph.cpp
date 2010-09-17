#include "config.h"
#include "ConstrainedSearch.h"
#include "ContigGraph.h"
#include "ContigPath.h"
#include "Estimate.h"
#include "Histogram.h"
#include "Iterator.h"
#include "Uncompress.h"
#include <algorithm> // for min
#include <cerrno>
#include <climits> // for UINT_MAX
#include <cmath>
#include <cstring> // for strerror
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <iterator>
#include <map>
#include <pthread.h>
#include <set>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

#define PROGRAM "SimpleGraph"

static const char VERSION_MESSAGE[] =
PROGRAM " (" PACKAGE_NAME ") " VERSION "\n"
"Written by Jared Simpson and Shaun Jackman.\n"
"\n"
"Copyright 2010 Canada's Michael Smith Genome Science Centre\n";

static const char USAGE_MESSAGE[] =
"Usage: " PROGRAM " [OPTION]... ADJ DIST\n"
"Find paths through contigs using distance estimates.\n"
"  ADJ   adjacency of the contigs\n"
"  DIST  distance estimates between the contigs\n"
"\n"
"  -k, --kmer=KMER_SIZE  k-mer size\n"
"      --max-cost=COST   maximum computational cost\n"
"  -o, --out=FILE        write result to FILE\n"
"  -j, --threads=THREADS use THREADS parallel threads [1]\n"
"      --scaffold        join contigs with Ns [default]\n"
"      --no-scaffold     do not scaffold\n"
"  -v, --verbose         display verbose output\n"
"      --help            display this help and exit\n"
"      --version         output version information and exit\n"
"\n"
"Report bugs to <" PACKAGE_BUGREPORT ">.\n";

namespace opt {
	unsigned k; // used by ContigProperties
	static unsigned threads = 1;
	static int scaffold = 1;
	static int verbose;
	static string out;
	int dot; // used by Estimate
}

static const char shortopts[] = "j:k:o:v";

enum { OPT_HELP = 1, OPT_VERSION, OPT_MAX_COST };

static const struct option longopts[] = {
	{ "kmer",        required_argument, NULL, 'k' },
	{ "max-cost",    required_argument, NULL, OPT_MAX_COST },
	{ "out",         required_argument, NULL, 'o' },
	{ "scaffold",    no_argument,       &opt::scaffold, 1 },
	{ "no-scaffold", no_argument,       &opt::scaffold, 0 },
	{ "threads",     required_argument,	NULL, 'j' },
	{ "verbose",     no_argument,       NULL, 'v' },
	{ "help",        no_argument,       NULL, OPT_HELP },
	{ "version",     no_argument,       NULL, OPT_VERSION },
	{ NULL, 0, NULL, 0 }
};

static void generatePathsThroughEstimates(const Graph& g,
		const string& estPath);

static void assert_open(ifstream& f, const string& p)
{
	if (f.is_open())
		return;
	cerr << p << ": " << strerror(errno) << endl;
	exit(EXIT_FAILURE);
}

int main(int argc, char** argv)
{
	bool die = false;
	for (int c; (c = getopt_long(argc, argv,
					shortopts, longopts, NULL)) != -1;) {
		istringstream arg(optarg != NULL ? optarg : "");
		switch (c) {
			case '?': die = true; break;
			case 'j': arg >> opt::threads; break;
			case 'k': arg >> opt::k; break;
			case OPT_MAX_COST: arg >> opt::maxCost; break;
			case 'o': arg >> opt::out; break;
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

	if (opt::out.empty()) {
		cerr << PROGRAM ": " << "missing -o,--out option\n";
		die = true;
	}

	if (argc - optind < 2) {
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

	string adjFile(argv[optind++]);
	string estFile(argv[optind++]);

	// Read the contig adjacency graph.
	ifstream fin(adjFile.c_str());
	assert_open(fin, adjFile);
	Graph g;
	fin >> g;
	assert(fin.eof());

	if (opt::verbose > 0) {
		unsigned v = g.num_vertices();
		unsigned e = g.num_edges();
		cout << "V=" << v << " E=" << e
			<< " E/V=" << (float)e / v << endl;

		// Print a histogram of the degree.
		Histogram h;
		typedef Graph::vertex_iterator vertex_iterator;
		std::pair<vertex_iterator, vertex_iterator>
			vit = vertices(g);
		for (vertex_iterator u = vit.first; u != vit.second; ++u)
			h.insert(out_degree(*u, g));
		cout <<
			"Degree: " << h.barplot() << "\n"
			"        01234" << endl;
	}

	// try to find paths that match the distance estimates
	generatePathsThroughEstimates(g, estFile);
}

static ostream& printConstraints(ostream& out, Constraints s)
{
	for (Constraints::const_iterator it = s.begin();
			it != s.end(); ++it)
		out << ' ' << it->first << ',' << it->second;
	return out;
}

/** Return the set of contigs that appear more than once in a single
 * solution.
 */
static set<ContigID> findRepeats(ContigID seed,
	const ContigPaths& solutions)
{
	set<ContigID> repeats;
	for (ContigPaths::const_iterator solIt = solutions.begin();
			solIt != solutions.end(); ++solIt) {
		map<ContigID, unsigned> count;
		count[seed]++;
		for (ContigPath::const_iterator it = solIt->begin();
				it != solIt->end(); ++it)
			count[ContigID(*it)]++;
		for (map<ContigID, unsigned>::const_iterator
				it = count.begin(); it != count.end(); ++it)
			if (it->second > 1)
				repeats.insert(it->first);
	}
	return repeats;
}

/** The fewest number of pairs in a distance estimate. */
static unsigned g_minNumPairs = UINT_MAX;

/** The fewest number of pairs used in a path. */
static unsigned g_minNumPairsUsed = UINT_MAX;

static struct {
	unsigned totalAttempted;
	unsigned uniqueEnd;
	unsigned noPossiblePaths;
	unsigned noValidPaths;
	unsigned repeat;
	unsigned multiEnd;
	unsigned tooManySolutions;
	unsigned tooComplex;
} stats;

/** Return the length of the specified path. */
static size_t calculatePathLength(const Graph& g,
		const ContigPath& path,
		size_t first = 0, size_t last = 0)
{
	size_t len = 0;
	if (first + last < path.size()) {
		for (ContigPath::const_iterator iter = path.begin() + first;
				iter != path.end() - last; ++iter)
			len += g[*iter].length;
	}
	assert(len > 0);
	return len;
}

/** Return an ambiguous path that agrees with all the given paths. */
static ContigPath constructAmbiguousPath(const Graph &g,
		const ContigPaths& solutions)
{
	typedef vector<ContigPath> ContigPaths;
	const ContigPaths& paths = solutions;

	// Find the size of the smallest path.
	const ContigPath& firstSol = paths.front();
	size_t min_len = firstSol.size();
	for (ContigPaths::const_iterator it = paths.begin() + 1;
			it != paths.end(); ++it)
		min_len = min(min_len, it->size());

	// Find the longest prefix.
	ContigPath vppath;
	size_t longestPrefix;
	bool commonPrefix = true;
	for (longestPrefix = 0;
			longestPrefix < min_len; longestPrefix++) {
		const ContigNode& common_path_node = firstSol[longestPrefix];
		for (ContigPaths::const_iterator solIter = paths.begin();
				solIter != paths.end(); ++solIter) {
			const ContigNode& pathnode = (*solIter)[longestPrefix];
			if (pathnode != common_path_node) {
				// Found the longest prefix.
				commonPrefix = false;
				break;
			}
		}
		if (!commonPrefix)
			break;
		vppath.push_back(common_path_node);
	}

	// Find the longest suffix.
	ContigPath vspath;
	size_t longestSuffix;
	bool commonSuffix = true;
	for (longestSuffix = 0;
			longestSuffix < min_len-longestPrefix; longestSuffix++) {
		const ContigNode& common_path_node
			= firstSol[firstSol.size()-longestSuffix-1];
		for (ContigPaths::const_iterator solIter = paths.begin();
				solIter != paths.end(); ++solIter) {
			const ContigNode& pathnode
				= (*solIter)[solIter->size()-longestSuffix-1];
			if (pathnode != common_path_node) {
				// Found the longest suffix.
				commonSuffix = false;
				break;
			}
		}
		if (!commonSuffix)
			break;
		vspath.push_back(common_path_node);
	}

	// Calculate the length of the longest path.
	unsigned maxLen = 0;
	ContigPaths::const_iterator longest = paths.end();
	for (ContigPaths::const_iterator it = paths.begin();
			it != paths.end(); ++it) {
		unsigned len = calculatePathLength(g, *it);
		if (len > maxLen) {
			maxLen = len;
			longest = it;
		}
	}
	assert(maxLen > 0);
	assert(longest != paths.end());

	unsigned numN = calculatePathLength(g, *longest,
			longestPrefix, longestSuffix);
	ContigPath out;
	out.reserve(vppath.size() + 1 + vspath.size());
	out.insert(out.end(), vppath.begin(), vppath.end());
	if (longestSuffix > 0) {
		out.push_back(ContigNode(numN, 'N'));
		out.insert(out.end(), vspath.rbegin(), vspath.rend());
	}
	return out;
}

/** Return a map of contig IDs to their distance along this path.
 * Repeat contigs, which would have more than one position, are not
 * represented in this map.
 */
void makeDistanceMap(const Graph& g,
		const ContigPath& path,
		map<ContigNode, int>& distances)
{
	size_t distance = 0;
	for (ContigPath::const_iterator it = path.begin();
			it != path.end(); ++it) {
		bool inserted = distances.insert(
				make_pair(*it, distance)).second;
		if (!inserted) {
			// Mark this contig as a repeat.
			distances[*it] = INT_MIN;
		}
		distance += g[*it].length;
	}

	// Remove the repeats.
	for (map<ContigNode, int>::iterator it = distances.begin();
			it != distances.end();)
		if (it->second == INT_MIN)
			distances.erase(it++);
		else
			++it;
}

/** Find a path for the specified distance estimates.
 * @param out [out] the solution path
 */
static void handleEstimate(const Graph& g,
		const EstimateRecord& er, bool dirIdx,
		ContigPath& out)
{
	if (er.estimates[dirIdx].empty())
		return;

	ostringstream vout_ss;
	ostream bitBucket(NULL);
	ostream& vout = opt::verbose > 0 ? vout_ss : bitBucket;
	vout << "\n* " << ContigNode(er.refID, dirIdx) << '\n';

	unsigned minNumPairs = UINT_MAX;
	// generate the reachable set
	Constraints constraints;
	for (EstimateVector::const_iterator iter
				= er.estimates[dirIdx].begin();
			iter != er.estimates[dirIdx].end(); ++iter) {
		minNumPairs = min(minNumPairs, iter->numPairs);

		// Translate the distances produced by the esimator into the
		// coordinate space used by the graph (a translation of k-1).
		constraints.push_back(Constraint(iter->contig,
					iter->distance + opt::k - 1
					+ allowedError(iter->stdDev)));
	}

	vout << "Constraints:";
	printConstraints(vout, constraints) << '\n';

	ContigPaths solutions;
	unsigned numVisited = 0;
	constrainedSearch(g, ContigNode(er.refID, dirIdx),
			constraints, solutions, numVisited);
	bool tooComplex = numVisited >= opt::maxCost;
	bool tooManySolutions = solutions.size() > opt::maxPaths;

	set<ContigID> repeats = findRepeats(er.refID, solutions);
	if (!repeats.empty()) {
		vout << "Repeats:";
		copy(repeats.begin(), repeats.end(),
				affix_ostream_iterator<ContigID>(vout, " "));
		vout << '\n';
	}

	unsigned numPossiblePaths = solutions.size();
	if (numPossiblePaths > 0)
		vout << "Paths: " << numPossiblePaths << '\n';

	for (ContigPaths::iterator solIter = solutions.begin();
			solIter != solutions.end();) {
		vout << *solIter << '\n';

		// Calculate the path distance to each node and see if
		// it is within the estimated distance.
		map<ContigNode, int> distanceMap;
		makeDistanceMap(g, *solIter, distanceMap);

		// Remove solutions whose distance estimates are not correct.
		unsigned validCount = 0, invalidCount = 0, ignoredCount = 0;
		for (EstimateVector::const_iterator iter
					= er.estimates[dirIdx].begin();
				iter != er.estimates[dirIdx].end(); ++iter) {
			vout << *iter << '\t';

			map<ContigNode, int>::iterator dmIter
				= distanceMap.find(iter->contig);
			if (dmIter == distanceMap.end()) {
				// This contig is a repeat.
				ignoredCount++;
				vout << "ignored\n";
				continue;
			}

			// translate distance by -overlap to match
			// coordinate space used by the estimate
			int actualDistance = dmIter->second - opt::k + 1;
			int diff = actualDistance - iter->distance;
			unsigned buffer = allowedError(iter->stdDev);
			bool invalid = (unsigned)abs(diff) > buffer;
			bool repeat = repeats.count(ContigID(iter->contig)) > 0;
			bool ignored = invalid && repeat;
			if (ignored)
				ignoredCount++;
			else if (invalid)
				invalidCount++;
			else
				validCount++;
			vout << "dist: " << actualDistance
				<< " diff: " << diff
				<< " buffer: " << buffer
				<< " n: " << iter->numPairs
				<< (ignored ? " ignored" : invalid ? " invalid" : "")
				<< '\n';
		}

		if (invalidCount == 0 && validCount > 0)
			++solIter;
		else
			solIter = solutions.erase(solIter);
	}

	vout << "Solutions: " << solutions.size();
	if (tooComplex)
		vout << " (too complex)";
	if (tooManySolutions)
		vout << " (too many solutions)";
	vout << '\n';

	ContigPaths::iterator bestSol = solutions.end();
	int minDiff = 999999;
	for (ContigPaths::iterator solIter = solutions.begin();
			solIter != solutions.end(); ++solIter) {
		map<ContigNode, int> distanceMap;
		makeDistanceMap(g, *solIter, distanceMap);
		int sumDiff = 0;
		for (EstimateVector::const_iterator iter
					= er.estimates[dirIdx].begin();
				iter != er.estimates[dirIdx].end(); ++iter) {
			if (repeats.count(ContigID(iter->contig)) > 0)
				continue;
			map<ContigNode, int>::iterator dmIter
				= distanceMap.find(iter->contig);
			assert(dmIter != distanceMap.end());
			int actualDistance = dmIter->second - opt::k + 1;
			int diff = actualDistance - iter->distance;
			sumDiff += abs(diff);
		}

		if (sumDiff < minDiff) {
			minDiff = sumDiff;
			bestSol = solIter;
		}

		size_t len = calculatePathLength(g, *solIter);
		vout << *solIter
			<< " length: " << len
			<< " sumdiff: " << sumDiff << '\n';
	}

	/** Lock the debugging stream. */
	static pthread_mutex_t coutMutex = PTHREAD_MUTEX_INITIALIZER;
	pthread_mutex_lock(&coutMutex);
	stats.totalAttempted++;
	g_minNumPairs = min(g_minNumPairs, minNumPairs);

	if (tooComplex) {
		stats.tooComplex++;
	} else if (tooManySolutions) {
		stats.tooManySolutions++;
	} else if (numPossiblePaths == 0) {
		stats.noPossiblePaths++;
	} else if (solutions.empty()) {
		stats.noValidPaths++;
	} else if (repeats.count(er.refID) > 0) {
		vout << "Repeat: " << ContigNode(er.refID, dirIdx) << '\n';
		stats.repeat++;
	} else if (solutions.size() > 1) {
		ContigPath path
			= constructAmbiguousPath(g, solutions);
		if (!path.empty()) {
			vout << path << '\n';
			if (opt::scaffold) {
				out.insert(out.end(), path.begin(), path.end());
				g_minNumPairsUsed
					= min(g_minNumPairsUsed, minNumPairs);
			}
		}
		stats.multiEnd++;
	} else {
		assert(solutions.size() == 1);
		assert(bestSol != solutions.end());
		out.insert(out.end(), bestSol->begin(), bestSol->end());
		stats.uniqueEnd++;
		g_minNumPairsUsed = min(g_minNumPairsUsed, minNumPairs);
	}
	cout << vout_ss.str();
	if (!out.empty())
		assert(!out.back().ambiguous());
	pthread_mutex_unlock(&coutMutex);
}

struct WorkerArg {
	istream* in;
	ostream* out;
	const Graph* graph;
	WorkerArg(istream* in, ostream* out, const Graph* g)
		: in(in), out(out), graph(g) { }
};

static void* worker(void* pArg)
{
	WorkerArg& arg = *static_cast<WorkerArg*>(pArg);
	for (;;) {
		/** Lock the input stream. */
		static pthread_mutex_t inMutex = PTHREAD_MUTEX_INITIALIZER;
		pthread_mutex_lock(&inMutex);
		EstimateRecord er;
		bool good = (*arg.in) >> er;
		pthread_mutex_unlock(&inMutex);
		if (!good)
			break;

		// Flip the anterior distance estimates.
		for (EstimateVector::iterator it = er.estimates[1].begin();
				it != er.estimates[1].end(); ++it)
			it->contig.flip();

		ContigPath path;
		handleEstimate(*arg.graph, er, true, path);
		path.reverseComplement();
		path.push_back(ContigNode(er.refID, false));
		handleEstimate(*arg.graph, er, false, path);
		if (path.size() > 1) {
			/** Lock the output stream. */
			static pthread_mutex_t outMutex
				= PTHREAD_MUTEX_INITIALIZER;
			pthread_mutex_lock(&outMutex);
			*arg.out << er.refID << '\t' << path << '\n';
			assert(arg.out->good());
			pthread_mutex_unlock(&outMutex);
		}
	}
	return NULL;
}

static void generatePathsThroughEstimates(const Graph& g,
		const string& estPath)
{
	ifstream inStream(estPath.c_str());
	assert_open(inStream, estPath);

	ofstream outStream(opt::out.c_str());
	assert(outStream.is_open());

	// Create the worker threads.
	vector<pthread_t> threads;
	threads.reserve(opt::threads);
	WorkerArg arg(&inStream, &outStream, &g);
	for (unsigned i = 0; i < opt::threads; i++) {
		pthread_t thread;
		pthread_create(&thread, NULL, worker, &arg);
		threads.push_back(thread);
	}

	// Wait for the worker threads to finish.
	for (vector<pthread_t>::const_iterator it = threads.begin();
			it != threads.end(); ++it) {
		void* status;
		pthread_join(*it, &status);
	}
	if (opt::verbose > 0)
		cout << '\n';

	cout <<
		"Total paths attempted: " << stats.totalAttempted << "\n"
		"Unique path: " << stats.uniqueEnd << "\n"
		"No possible paths: " << stats.noPossiblePaths << "\n"
		"No valid paths: " << stats.noValidPaths << "\n"
		"Repetitive: " << stats.repeat << "\n"
		"Multiple valid paths: " << stats.multiEnd << "\n"
		"Too many solutions: " << stats.tooManySolutions << "\n"
		"Too complex: " << stats.tooComplex << "\n";

	inStream.close();
	outStream.close();

	cout << "\n"
		"The minimum number of pairs in a distance estimate is "
		<< g_minNumPairs << ".\n";
	if (g_minNumPairsUsed != UINT_MAX) {
		cout << "The minimum number of pairs used in a path is "
			<< g_minNumPairsUsed << ".\n";
		if (g_minNumPairs < g_minNumPairsUsed)
			cout << "Consider increasing the number of pairs "
				"threshold paramter, n, to " << g_minNumPairsUsed
				<< ".\n";
	}
}
