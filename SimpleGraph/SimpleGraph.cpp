#include "config.h"
#include "ContigGraph.h"
#include "ContigPath.h"
#include "DirectedGraphImpl.h"
#include "PairUtils.h"
#include "Uncompress.h"
#include <cmath>
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <iterator>
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


// Typedefs
typedef std::vector<LinearNumKey> LinearNumVector;
typedef std::vector<LinearNumVector > LinearNumKeyVector;

// FUNCTORS
struct SimpleDataCost
{
	SimpleDataCost(size_t kmer) : m_overlap(kmer - 1) { }
	
	size_t cost(const SimpleContigData& data) { return data.length - m_overlap; } 
	size_t m_overlap;
};


// FUNCTIONS
void generatePathsThroughEstimates(SimpleContigGraph* pContigGraph, std::string estFileName, int kmer, int maxCost);
void constructContigPath(const SimpleContigGraph::VertexPath& vertexPath, ContigPath& contigPath);
void outputContigPath(std::ofstream& outStream, LinearNumKey refNode, extDirection dir, const ContigPath& contigPath);

//
//
//
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

	std::cout << "Adj file: " << adjFile
		<< " Estimate File: " << estFile
		<< " len file: " << lenFile
		<< " kmer " << opt::k
		<< " max cost " << opt::maxCost
		<< endl;

	// Load the graph from the adjacency file
	SimpleContigGraph contigGraph;
	loadGraphFromAdjFile(&contigGraph, lenFile, adjFile);

	// try to find paths that match the distance estimates
	generatePathsThroughEstimates(&contigGraph, estFile,
			opt::k, opt::maxCost);
}

template<typename K, typename D>
ostream& printMap(const map<K,D>& s)
{
	cout << "{ ";
	for (typename map<K,D>::const_iterator iter = s.begin();
			iter != s.end(); ++iter)
		cout << iter->first << ',' << iter->second << ' ';
	return cout << '}';
}

template<typename K>
ostream& printPath(const vector<K>& s)
{
	cout << "[ ";
	copy(s.begin(), s.end(),
			ostream_iterator<K>(cout, " "));
	return cout << ']';
}

void generatePathsThroughEstimates(SimpleContigGraph* pContigGraph, std::string estFileName, int kmer, int maxCost)
{
	int totalAttempted = 0;
	int noPossiblePaths = 0;
	int nopathEnd = 0;
	int uniqueEnd = 0;
	int multiEnd = 0;

	(void)pContigGraph;
	std::ifstream inStream(estFileName.c_str());
	std::ofstream outStream(opt::out.c_str());
	assert(outStream.is_open());
	SimpleDataCost costFunctor(kmer);

	for (EstimateRecord er; inStream >> er;) {
		for(size_t dirIdx = 0; dirIdx <= 1; ++dirIdx)
		{
			bool gDebugPrint = opt::verbose > 0;
		  if(gDebugPrint) std::cout << "****Processing " << er.refID << " dir: " << dirIdx << "****\n";
			// generate the reachable set
			SimpleContigGraph::KeyConstraintMap constraintMap;
			
			for(EstimateVector::iterator iter = er.estimates[dirIdx].begin(); iter != er.estimates[dirIdx].end(); ++iter)
			{
				// Translate the distances produced by the esimator into the coordinate space
				// used by the graph (a translation of m_overlap)
				int translatedDistance = iter->distance + costFunctor.m_overlap;
				unsigned distanceBuffer = allowedError(iter->stdDev);
				
				Constraint nc;
				nc.distance = translatedDistance  + distanceBuffer;
				nc.isRC = iter->isRC;

				if(gDebugPrint) std::cout << "Adding Constraint " << iter->nID << ' ' << nc << "\n";

				constraintMap[iter->nID] = nc;
			}

			if(gDebugPrint)
			{
				std::cout << "Constraints: ";
				printMap(constraintMap);
				std::cout << "\n";
			}

			// Generate the paths
			SimpleContigGraph::FeasiblePaths solutions;
			
			// The value at which to abort trying to find a path
			const int maxNumPaths = 200;
			int numVisited = 0;
			
			// The number of nodes to explore before bailing out
			pContigGraph->findSuperpaths(er.refID, (extDirection)dirIdx, constraintMap, solutions, costFunctor, maxNumPaths, maxCost, numVisited);
			
			// Select the best path
			//std::cout << "Computational cost " << numVisited << "\n";
			if (gDebugPrint)
				std::cout << "Possible solutions: "
					<< solutions.size() << "\n";

			totalAttempted++;
			if (solutions.empty()) {
				noPossiblePaths++;
				continue;
			}

			for (SimpleContigGraph::FeasiblePaths::iterator
					solIter = solutions.begin();
					solIter != solutions.end();) {
				if (gDebugPrint) {
					printPath(*solIter);
					std::cout << '\n';
				}

				// Calculate the path distance to each node and see if
				// it is within the estimated distance.
				std::map<LinearNumKey, int> distanceMap;
				pContigGraph->makeDistanceMap(*solIter,
						costFunctor, distanceMap);

				// Filter out potentially bad hits
				bool validPath = true;
				for (EstimateVector::iterator
						iter = er.estimates[dirIdx].begin();
						iter != er.estimates[dirIdx].end(); ++iter) {
					if (gDebugPrint)
						std::cout << "Estimate " << *iter << "\n";

					std::map<LinearNumKey, int>::iterator
						dmIter = distanceMap.find(iter->nID);
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
						std::cout
							<< "Actual dist: " << actualDistance
							<< " diff: " << diff
							<< " buffer: " << buffer
							<< " n: " << iter->numPairs
							<< (invalid ? " invalid" : "")
							<< "\n";
				}

				if (validPath)
					++solIter;
				else
					solIter = solutions.erase(solIter);
			}

			if (gDebugPrint)
				std::cout << "Solutions: "
					<< solutions.size() << "\n";

			SimpleContigGraph::FeasiblePaths::iterator bestSol
				= solutions.end();
			int minDiff = 999999;
			for (SimpleContigGraph::FeasiblePaths::iterator
					solIter = solutions.begin();
					solIter != solutions.end(); ++solIter) {
				size_t len = pContigGraph->calculatePathLength(
						*solIter, costFunctor);

				if (gDebugPrint) {
					printPath(*solIter);
					std::cout << " length: " << len;
				}

				std::map<LinearNumKey, int> distanceMap;
				pContigGraph->makeDistanceMap(*solIter,
						costFunctor, distanceMap);
				int sumDiff = 0;
				for (EstimateVector::iterator iter
						= er.estimates[dirIdx].begin();
						iter != er.estimates[dirIdx].end(); ++iter) {
					std::map<LinearNumKey, int>::iterator dmIter
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
					std::cout << " sumdiff: " << sumDiff << "\n";
			}

			if (solutions.empty()) {
				nopathEnd++;
			} else if (solutions.size() > 1) {
				multiEnd++;
			} else {
				assert(solutions.size() == 1);
				assert(bestSol != solutions.end());
				uniqueEnd++;
				ContigPath cPath;
				constructContigPath(*bestSol, cPath);
				outputContigPath(outStream,
						er.refID, (extDirection)dirIdx, cPath);
				assert(outStream.good());
			}
		}
	}
	
	cout <<
		"Total paths attempted: " << totalAttempted << "\n"
		"No possible paths: " << noPossiblePaths << "\n"
		"No valid paths: " << nopathEnd << "\n"
		"Multiple valid paths: " << multiEnd << "\n"
		"Unique path: " << uniqueEnd << "\n";
	
	inStream.close();
	outStream.close();
}

//
// Convert the vertext path to a contig path
// The difference between the two is that the complementness of a vertex path is with respect to the previous node in the path
// but with a contig path it is with respect to the root contig. Also the contigPath uses true contig ids
// 
void constructContigPath(const SimpleContigGraph::VertexPath& vertexPath, ContigPath& contigPath)
{
	bool flip = false;
	for(SimpleContigGraph::VertexPath::const_iterator iter = vertexPath.begin(); iter != vertexPath.end(); ++iter)
	{
		flip = flip ^ iter->isRC;
		MergeNode mn = {iter->key, flip};
		contigPath.appendNode(mn);
	}
}

void outputContigPath(std::ofstream& outStream, LinearNumKey refNode, extDirection dir, const ContigPath& contigPath)
{
  //std::cout << "Found path for " << refNode << "\n";
	outStream << "@ " << refNode << "," << dir << " -> ";
	outStream << contigPath << "\n";
}
