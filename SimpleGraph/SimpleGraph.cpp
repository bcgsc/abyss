#include "config.h"
#include "ContigPath.h"
#include "Estimate.h"
#include "IOUtil.h"
#include "Uncompress.h"
#include "Graph/ConstrainedSearch.h"
#include "Graph/ContigGraph.h"
#include "Graph/ContigGraphAlgorithms.h"
#include "Graph/GraphIO.h"
#include "Graph/GraphUtil.h"
#include <algorithm> // for min
#include <climits> // for UINT_MAX
#include <cmath>
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
"Copyright 2012 Canada's Michael Smith Genome Science Centre\n";

static const char USAGE_MESSAGE[] =
"Usage: " PROGRAM " [OPTION]... ADJ DIST\n"
"Find paths through contigs using distance estimates.\n"
"  ADJ   adjacency of the contigs\n"
"  DIST  distance estimates between the contigs\n"
"\n"
"  -k, --kmer=KMER_SIZE  k-mer size\n"
"  -d, --dist-error=N    acceptable error of a distance estimate\n"
"                        default is 6 bp\n"
"      --max-cost=COST   maximum computational cost\n"
"  -o, --out=FILE        write result to FILE\n"
"  -j, --threads=THREADS use THREADS parallel threads [1]\n"
"      --extend          extend unambiguous paths\n"
"      --no-extend       do not extend unambiguous paths [default]\n"
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
	static int extend;
	static int scaffold = 1;
	static int verbose;
	static string out;

	/** The acceptable error of a distance estimate. */
	unsigned distanceError = 6;

 	/** Output format */
 	int format = DIST; // used by Estimate
}

static const char shortopts[] = "d:j:k:o:v";

enum { OPT_HELP = 1, OPT_VERSION, OPT_MAX_COST };

static const struct option longopts[] = {
	{ "kmer",        required_argument, NULL, 'k' },
	{ "dist-error",  required_argument, NULL, 'd' },
	{ "max-cost",    required_argument, NULL, OPT_MAX_COST },
	{ "out",         required_argument, NULL, 'o' },
	{ "extend",      no_argument,       &opt::extend, 1 },
	{ "no-extend",   no_argument,       &opt::extend, 0 },
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

int main(int argc, char** argv)
{
	bool die = false;
	for (int c; (c = getopt_long(argc, argv,
					shortopts, longopts, NULL)) != -1;) {
		istringstream arg(optarg != NULL ? optarg : "");
		switch (c) {
			case '?': die = true; break;
			case 'd': arg >> opt::distanceError; break;
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
	if (opt::verbose > 0)
		cerr << "Reading `" << adjFile << "'..." << endl;
	ifstream fin(adjFile.c_str());
	assert_good(fin, adjFile);
	Graph g;
	fin >> g;
	assert(fin.eof());

	if (opt::verbose > 0)
		printGraphStats(cout, g);

	// try to find paths that match the distance estimates
	generatePathsThroughEstimates(g, estFile);
}

/** Print a set of constraints. */
static ostream& printConstraints(ostream& out,
		const Graph& g, const Constraints& s)
{
	for (Constraints::const_iterator it = s.begin();
			it != s.end(); ++it)
		out << ' ' << get(vertex_name, g, it->first)
			<< ',' << it->second;
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
			count[it->contigIndex()]++;
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

typedef graph_traits<Graph>::vertex_descriptor vertex_descriptor;

/** Return the distance from vertex u to v. */
static int getDistance(const Graph& g,
		vertex_descriptor u, vertex_descriptor v)
{
	typedef graph_traits<Graph>::edge_descriptor edge_descriptor;
	pair<edge_descriptor, bool> e = edge(u, v, g);
	assert(e.second);
	return g[e.first].distance;
}

/** Return the length of the specified path in k-mer. */
static unsigned calculatePathLength(const Graph& g,
		const ContigNode& origin,
		const ContigPath& path, size_t prefix = 0, size_t suffix = 0)
{
	if (prefix + suffix == path.size())
		return 0;
	assert(prefix + suffix < path.size());
	int length = addProp(g, path.begin() + prefix,
			path.end() - suffix).length;

	// Account for the overlap on the left.
	vertex_descriptor u = prefix == 0 ? origin : path[prefix - 1];
	length += getDistance(g, u, path[prefix]);
	assert(length > 0);
	return length;
}

/** Compare the lengths of two paths. */
struct ComparePathLength
		: binary_function<ContigPath, ContigPath, bool>
{
	ComparePathLength(const Graph& g, const ContigNode& origin)
		: m_g(g), m_origin(origin) { }
	bool operator()(const ContigPath& a, const ContigPath& b) const {
		unsigned lenA = calculatePathLength(m_g, m_origin, a);
		unsigned lenB = calculatePathLength(m_g, m_origin, b);
		return lenA < lenB
			|| lenA == lenB && a.size() < b.size();
	}
  private:
	const Graph& m_g;
	const ContigNode& m_origin;
};

/** Return an ambiguous path that agrees with all the given paths. */
static ContigPath constructAmbiguousPath(const Graph &g,
		const ContigNode& origin, const ContigPaths& paths)
{
	assert(!paths.empty());

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

	ContigPath out;
	out.reserve(vppath.size() + 1 + vspath.size());
	out.insert(out.end(), vppath.begin(), vppath.end());
	if (longestSuffix > 0) {
		const ContigPath& longestPath(
				*max_element(paths.begin(), paths.end(),
					ComparePathLength(g, origin)));
		unsigned length = calculatePathLength(g, origin, longestPath,
				longestPrefix, longestSuffix);

		// Account for the overlap on the right.
		int dist = length + getDistance(g,
				longestSuffix == longestPath.size() ? origin
				: *(longestPath.rbegin() + longestSuffix),
				*(longestPath.rbegin() + longestSuffix - 1));

		// Add k-1 because it is the convention.
		int numN = dist + opt::k - 1;
		assert(numN > 0);

		out.push_back(ContigNode(numN, 'N'));
		out.insert(out.end(), vspath.rbegin(), vspath.rend());
	}
	return out;
}

/** Return a map of contig IDs to their distance along this path.
 * Repeat contigs, which would have more than one position, are not
 * represented in this map.
 */
map<ContigNode, int> makeDistanceMap(const Graph& g,
		const ContigNode& origin, const ContigPath& path)
{
	map<ContigNode, int> distances;
	int distance = 0;
	for (ContigPath::const_iterator it = path.begin();
			it != path.end(); ++it) {
		vertex_descriptor u = it == path.begin() ? origin : *(it - 1);
		vertex_descriptor v = *it;
		distance += getDistance(g, u, v);

		bool inserted = distances.insert(
				make_pair(v, distance)).second;
		if (!inserted) {
			// Mark this contig as a repeat.
			distances[v] = INT_MIN;
		}

		distance += g[v].length;
	}

	// Remove the repeats.
	for (map<ContigNode, int>::iterator it = distances.begin();
			it != distances.end();)
		if (it->second == INT_MIN)
			distances.erase(it++);
		else
			++it;
	return distances;
}

/** Print a distance map. */
static void printDistanceMap(ostream& out, const Graph& g,
		const ContigNode& u, const ContigPath& path)
{
	typedef map<ContigNode, int> DistanceMap;
	DistanceMap distanceMap = makeDistanceMap(g, u, path);
	for (DistanceMap::const_iterator it = distanceMap.begin();
			it != distanceMap.end(); ++it)
		out << get(edge_name, g, make_pair(u, it->first))
			<< " [d=" << it->second << "]\n";
}

typedef std::vector<std::pair<ContigNode, DistanceEst> > Estimates;

/** Find a path for the specified distance estimates.
 * @param out [out] the solution path
 */
static void handleEstimate(const Graph& g,
		const EstimateRecord& er, bool dirIdx,
		ContigPath& out)
{
	if (er.estimates[dirIdx].empty())
		return;

	ContigNode origin(er.refID, dirIdx);
	ostringstream vout_ss;
	ostream bitBucket(NULL);
	ostream& vout = opt::verbose > 0 ? vout_ss : bitBucket;
	vout << "\n* " << get(vertex_name, g, origin) << '\n';

	unsigned minNumPairs = UINT_MAX;
	// generate the reachable set
	Constraints constraints;
	for (Estimates::const_iterator iter
				= er.estimates[dirIdx].begin();
			iter != er.estimates[dirIdx].end(); ++iter) {
		ContigNode v = iter->first;
		const DistanceEst& ep = iter->second;
		minNumPairs = min(minNumPairs, ep.numPairs);
		constraints.push_back(Constraint(v,
					ep.distance + allowedError(ep.stdDev)));
	}

	vout << "Constraints:";
	printConstraints(vout, g, constraints) << '\n';

	ContigPaths solutions;
	unsigned numVisited = 0;
	constrainedSearch(g, origin, constraints, solutions, numVisited);
	bool tooComplex = numVisited >= opt::maxCost;
	bool tooManySolutions = solutions.size() > opt::maxPaths;

	set<ContigID> repeats = findRepeats(er.refID, solutions);
	if (!repeats.empty()) {
		vout << "Repeats:";
		for (set<ContigID>::const_iterator it = repeats.begin();
				it != repeats.end(); ++it)
			vout << ' ' << get(g_contigNames, *it);
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
		map<ContigNode, int> distanceMap
			= makeDistanceMap(g, origin, *solIter);

		// Remove solutions whose distance estimates are not correct.
		unsigned validCount = 0, invalidCount = 0, ignoredCount = 0;
		for (Estimates::const_iterator iter
					= er.estimates[dirIdx].begin();
				iter != er.estimates[dirIdx].end(); ++iter) {
			ContigNode v = iter->first;
			const DistanceEst& ep = iter->second;
			vout << get(vertex_name, g, v) << ',' << ep << '\t';

			map<ContigNode, int>::iterator dmIter
				= distanceMap.find(v);
			if (dmIter == distanceMap.end()) {
				// This contig is a repeat.
				ignoredCount++;
				vout << "ignored\n";
				continue;
			}

			// translate distance by -overlap to match
			// coordinate space used by the estimate
			int actualDistance = dmIter->second;
			int diff = actualDistance - ep.distance;
			unsigned buffer = allowedError(ep.stdDev);
			bool invalid = (unsigned)abs(diff) > buffer;
			bool repeat = repeats.count(v.contigIndex()) > 0;
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
				<< " n: " << ep.numPairs
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
		map<ContigNode, int> distanceMap
			= makeDistanceMap(g, origin, *solIter);
		int sumDiff = 0;
		for (Estimates::const_iterator iter
					= er.estimates[dirIdx].begin();
				iter != er.estimates[dirIdx].end(); ++iter) {
			ContigNode v = iter->first;
			const DistanceEst& ep = iter->second;
			if (repeats.count(v.contigIndex()) > 0)
				continue;
			map<ContigNode, int>::iterator dmIter
				= distanceMap.find(v);
			assert(dmIter != distanceMap.end());
			int actualDistance = dmIter->second;
			int diff = actualDistance - ep.distance;
			sumDiff += abs(diff);
		}

		if (sumDiff < minDiff) {
			minDiff = sumDiff;
			bestSol = solIter;
		}

		vout << *solIter
			<< " length: " << calculatePathLength(g, origin, *solIter)
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
		vout << "Repeat: " << get(vertex_name, g, origin) << '\n';
		stats.repeat++;
	} else if (solutions.size() > 1) {
		ContigPath path
			= constructAmbiguousPath(g, origin, solutions);
		if (!path.empty()) {
			if (opt::extend)
				extend(g, path.back(), back_inserter(path));
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
		ContigPath& path = *bestSol;
		if (opt::verbose > 1)
			printDistanceMap(vout, g, origin, path);
		if (opt::extend)
			extend(g, path.back(), back_inserter(path));
		out.insert(out.end(), path.begin(), path.end());
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
		for (Estimates::iterator it = er.estimates[1].begin();
				it != er.estimates[1].end(); ++it)
			it->first ^= 1;

		ContigPath path;
		handleEstimate(*arg.graph, er, true, path);
		reverseComplement(path.begin(), path.end());
		path.push_back(ContigNode(er.refID, false));
		handleEstimate(*arg.graph, er, false, path);
		if (path.size() > 1) {
			/** Lock the output stream. */
			static pthread_mutex_t outMutex
				= PTHREAD_MUTEX_INITIALIZER;
			pthread_mutex_lock(&outMutex);
			*arg.out << get(g_contigNames, er.refID)
				<< '\t' << path << '\n';
			assert(arg.out->good());
			pthread_mutex_unlock(&outMutex);
		}
	}
	return NULL;
}

static void generatePathsThroughEstimates(const Graph& g,
		const string& estPath)
{
	if (opt::verbose > 0)
		cerr << "Reading `" << estPath << "'..." << endl;
	ifstream inStream(estPath.c_str());
	assert_good(inStream, estPath);

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
