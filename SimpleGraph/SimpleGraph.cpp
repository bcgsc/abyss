#include "config.h"
#include "ContigGraph.h"
#include "ContigPath.h"
#include "DirectedGraphImpl.h"
#include "PairUtils.h"
#include "Uncompress.h"
#include <algorithm> // for min
#include <climits> // for UINT_MAX
#include <cmath>
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <iterator>
#include <pthread.h>
#include <sstream>
#include <vector>

using namespace std;

#define PROGRAM "SimpleGraph"

static const char VERSION_MESSAGE[] =
PROGRAM " (" PACKAGE_NAME ") " VERSION "\n"
"Written by Jared Simpson and Shaun Jackman.\n"
"\n"
"Copyright 2009 Canada's Michael Smith Genome Science Centre\n";

static const char USAGE_MESSAGE[] =
"Usage: " PROGRAM " [OPTION]... ADJ DIST\n"
"Find paths through contigs using distance estimates.\n"
"  ADJ   adjacency of the contigs\n"
"  DIST  distance estimates between the contigs\n"
"\n"
"  -k, --kmer=KMER_SIZE  k-mer size\n"
"      --max-cost=COST   maximum computational cost\n"
"  -o, --out=FILE        write result to FILE\n"
"  -v, --verbose         display verbose output\n"
"      --help            display this help and exit\n"
"      --version         output version information and exit\n"
"\n"
"Report bugs to <" PACKAGE_BUGREPORT ">.\n";

namespace opt {
	static unsigned k;
	static unsigned maxCost = 100000;
	static unsigned threads = 1;
	static int verbose;
	static string out;
}

static const char shortopts[] = "k:o:v";

enum { OPT_HELP = 1, OPT_VERSION, OPT_MAX_COST };

static const struct option longopts[] = {
	{ "kmer",        required_argument, NULL, 'k' },
	{ "max-cost",    required_argument, NULL, OPT_MAX_COST },
	{ "out",         required_argument, NULL, 'o' },
	{ "verbose",     no_argument,       NULL, 'v' },
	{ "help",        no_argument,       NULL, OPT_HELP },
	{ "version",     no_argument,       NULL, OPT_VERSION },
	{ NULL, 0, NULL, 0 }
};

typedef vector<LinearNumKey> LinearNumVector;
typedef vector<LinearNumVector > LinearNumKeyVector;

struct SimpleDataCost
{
	SimpleDataCost(size_t kmer) : m_overlap(kmer - 1) { }
	
	size_t cost(const SimpleContigData& data) { return data.length - m_overlap; } 
	size_t m_overlap;
};

static void generatePathsThroughEstimates(
		SimpleContigGraph* pContigGraph, string estFileName);
static void constructContigPath(
		const SimpleContigGraph::VertexPath& vertexPath,
		ContigPath& contigPath);

int main(int argc, char** argv)
{
	bool die = false;
	for (int c; (c = getopt_long(argc, argv,
					shortopts, longopts, NULL)) != -1;) {
		istringstream arg(optarg != NULL ? optarg : "");
		switch (c) {
			case '?': die = true; break;
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
	const string& lenFile = adjFile;

	// Load the graph from the adjacency file
	SimpleContigGraph contigGraph;
	loadGraphFromAdjFile(&contigGraph, lenFile, adjFile);
	if (opt::verbose > 0)
		cerr << "Vertices: " << contigGraph.getNumVertices()
			<< " Edges: " << contigGraph.countEdges() << endl;

	// try to find paths that match the distance estimates
	generatePathsThroughEstimates(&contigGraph, estFile);
}

template<typename K, typename D>
ostream& printMap(const map<K,D>& s)
{
	for (typename map<K,D>::const_iterator iter = s.begin();
			iter != s.end(); ++iter)
		cout << ' ' << iter->first << iter->second;
	return cout;
}

template<typename K>
ostream& printPath(const vector<K>& s)
{
	assert(!s.empty());
	copy(s.begin(), s.end()-1,
			ostream_iterator<K>(cout, " "));
	return cout << s.back();
}

static struct {
	unsigned totalAttempted;
	unsigned noPossiblePaths;
	unsigned nopathEnd;
	unsigned uniqueEnd;
	unsigned multiEnd;
} stats;

/** The fewest number of pairs in a distance estimate. */
static unsigned minNumPairs = UINT_MAX;

/** The fewest number of pairs used in a path. */
static unsigned minNumPairsUsed = UINT_MAX;

static void handleEstimate(const EstimateRecord& er, unsigned dirIdx,
		const SimpleContigGraph* pContigGraph, ostream& outStream)
{
	bool gDebugPrint = opt::verbose > 0;
	if (gDebugPrint)
		cout << "\n"
			"* " << er.refID << (dirIdx ? '-' : '+') << '\n';

	SimpleDataCost costFunctor(opt::k);
	// generate the reachable set
	SimpleContigGraph::KeyConstraintMap constraintMap;
	for (EstimateVector::const_iterator iter
				= er.estimates[dirIdx].begin();
			iter != er.estimates[dirIdx].end(); ++iter) {
		// Translate the distances produced by the esimator into the
		// coordinate space used by the graph (a translation of
		// m_overlap)
		int translatedDistance = iter->distance
			+ costFunctor.m_overlap;
		unsigned distanceBuffer = allowedError(iter->stdDev);

		Constraint nc;
		nc.distance = translatedDistance  + distanceBuffer;
		nc.isRC = iter->isRC;
		constraintMap[iter->nID] = nc;

		minNumPairs = min(minNumPairs, iter->numPairs);
	}

	if (gDebugPrint) {
		cout << "Constraints:";
		printMap(constraintMap);
		cout << '\n';
	}

	// Generate the paths
	SimpleContigGraph::FeasiblePaths solutions;

	// The value at which to abort trying to find a path
	const int maxNumPaths = 200;
	int numVisited = 0;

	// The number of nodes to explore before bailing out
	pContigGraph->findSuperpaths(er.refID, (extDirection)dirIdx,
			constraintMap, solutions, costFunctor, maxNumPaths,
			opt::maxCost, numVisited);

	if (gDebugPrint)
		cout << "Paths: " << solutions.size() << '\n';

	/** Lock the global variable stats. */
	static pthread_mutex_t statsMutex = PTHREAD_MUTEX_INITIALIZER;

	pthread_mutex_lock(&statsMutex);
	stats.totalAttempted++;
	pthread_mutex_unlock(&statsMutex);

	if (solutions.empty()) {
		pthread_mutex_lock(&statsMutex);
		stats.noPossiblePaths++;
		pthread_mutex_unlock(&statsMutex);
		return;
	}

	for (SimpleContigGraph::FeasiblePaths::iterator solIter
				= solutions.begin(); solIter != solutions.end();) {
		if (gDebugPrint) {
			printPath(*solIter);
			cout << '\n';
		}

		// Calculate the path distance to each node and see if
		// it is within the estimated distance.
		map<LinearNumKey, int> distanceMap;
		pContigGraph->makeDistanceMap(*solIter,
				costFunctor, distanceMap);

		// Filter out potentially bad hits
		bool validPath = true;
		for (EstimateVector::const_iterator iter
					= er.estimates[dirIdx].begin();
				iter != er.estimates[dirIdx].end(); ++iter) {
			if (gDebugPrint)
				cout << *iter << '\t';

			map<LinearNumKey, int>::iterator dmIter
				= distanceMap.find(iter->nID);
			assert(dmIter != distanceMap.end());
			// translate distance by -overlap to match
			// coordinate space used by the estimate
			int actualDistance
				= dmIter->second - costFunctor.m_overlap;
			int diff = actualDistance - iter->distance;
			unsigned buffer = allowedError(iter->stdDev);
			bool invalid = (unsigned)abs(diff) > buffer;
			if (invalid)
				validPath = false;
			if (gDebugPrint)
				cout << "dist: " << actualDistance
					<< " diff: " << diff
					<< " buffer: " << buffer
					<< " n: " << iter->numPairs
					<< (invalid ? " invalid" : "")
					<< '\n';
		}

		if (validPath)
			++solIter;
		else
			solIter = solutions.erase(solIter);
	}

	if (gDebugPrint)
		cout << "Solutions: " << solutions.size() << '\n';

	SimpleContigGraph::FeasiblePaths::iterator bestSol
		= solutions.end();
	int minDiff = 999999;
	for (SimpleContigGraph::FeasiblePaths::iterator solIter
				= solutions.begin();
			solIter != solutions.end(); ++solIter) {
		size_t len = pContigGraph->calculatePathLength(
				*solIter, costFunctor);

		if (gDebugPrint) {
			printPath(*solIter);
			cout << " length: " << len;
		}

		map<LinearNumKey, int> distanceMap;
		pContigGraph->makeDistanceMap(*solIter,
				costFunctor, distanceMap);
		int sumDiff = 0;
		for (EstimateVector::const_iterator iter
					= er.estimates[dirIdx].begin();
				iter != er.estimates[dirIdx].end(); ++iter) {
			map<LinearNumKey, int>::iterator dmIter
				= distanceMap.find(iter->nID);
			if (dmIter != distanceMap.end()) {
				int actualDistance = dmIter->second
					- costFunctor.m_overlap;
				int diff = actualDistance - iter->distance;
				sumDiff += abs(diff);
			}
		}

		if (sumDiff < minDiff) {
			minDiff = sumDiff;
			bestSol = solIter;
		}

		if (gDebugPrint)
			cout << " sumdiff: " << sumDiff << '\n';
	}

	if (solutions.empty()) {
		pthread_mutex_lock(&statsMutex);
		stats.nopathEnd++;
		pthread_mutex_unlock(&statsMutex);
	} else if (solutions.size() > 1) {
		pthread_mutex_lock(&statsMutex);
		stats.multiEnd++;
		pthread_mutex_unlock(&statsMutex);
	} else {
		assert(solutions.size() == 1);
		assert(bestSol != solutions.end());
		pthread_mutex_lock(&statsMutex);
		stats.uniqueEnd++;
		pthread_mutex_unlock(&statsMutex);
		ContigPath cPath;
		constructContigPath(*bestSol, cPath);

		/** Lock the output stream. */
		static pthread_mutex_t outMutex = PTHREAD_MUTEX_INITIALIZER;
		pthread_mutex_lock(&outMutex);
		outStream << "@ " << g_contigIDs.key(er.refID) << ','
			<< dirIdx << " -> " << cPath << '\n';
		assert(outStream.good());
		pthread_mutex_unlock(&outMutex);
	}

	const EstimateVector& v = er.estimates[dirIdx];
	for (EstimateVector::const_iterator it = v.begin();
			it != v.end(); ++it)
		minNumPairsUsed = min(minNumPairsUsed, it->numPairs);
}

struct WorkerArg {
	istream* in;
	ostream* out;
	const SimpleContigGraph* graph;
	WorkerArg(istream* in, ostream* out, const SimpleContigGraph* g)
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
		for (unsigned dirIdx = 0; dirIdx <= 1; ++dirIdx)
			handleEstimate(er, dirIdx, arg.graph, *arg.out);
	}
	return NULL;
}

static void generatePathsThroughEstimates(
		SimpleContigGraph* pContigGraph, string estFileName)
{
	ifstream inStream(estFileName.c_str());
	ofstream outStream(opt::out.c_str());
	assert(outStream.is_open());

	// Create the worker threads.
	vector<pthread_t> threads;
	threads.reserve(opt::threads);
	WorkerArg arg(&inStream, &outStream, pContigGraph);
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

	cout <<
		"Total paths attempted: " << stats.totalAttempted << "\n"
		"No possible paths: " << stats.noPossiblePaths << "\n"
		"No valid paths: " << stats.nopathEnd << "\n"
		"Multiple valid paths: " << stats.multiEnd << "\n"
		"Unique path: " << stats.uniqueEnd << "\n";

	inStream.close();
	outStream.close();

	cout << "\n"
		"The minimum number of pairs in a distance estimate is "
		<< minNumPairs << ".\n";
	if (minNumPairsUsed != UINT_MAX) {
		cout << "The minimum number of pairs used in a path is "
			<< minNumPairsUsed << ".\n";
		if (minNumPairs < minNumPairsUsed)
			cout << "Consider increasing the number of pairs "
				"threshold paramter, n, to " << minNumPairsUsed
				<< ".\n";
	}
}

//
// Convert the vertext path to a contig path
// The difference between the two is that the complementness of a vertex path is with respect to the previous node in the path
// but with a contig path it is with respect to the root contig. Also the contigPath uses true contig ids
//
static void constructContigPath(
		const SimpleContigGraph::VertexPath& vertexPath,
		ContigPath& contigPath)
{
	bool flip = false;
	for(SimpleContigGraph::VertexPath::const_iterator iter = vertexPath.begin(); iter != vertexPath.end(); ++iter)
	{
		flip = flip ^ iter->isRC;
		MergeNode mn = {iter->key, flip};
		contigPath.appendNode(mn);
	}
}
