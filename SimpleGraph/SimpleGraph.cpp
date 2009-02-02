#include <cmath>
#include <iostream>
#include <sstream>
#include <vector>
#include "SequenceCollectionHash.h"
#include "AssemblyAlgorithms.h"
#include "Options.h"
#include "FastaReader.h"
#include "Stats.h"
#include "DirectedGraph.h"
#include "PairUtils.h"
#include "ContigPath.h"
#include "SetOperations.h"

using namespace std;

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

static bool gDebugPrint = false;

//
//
//
int main(int argc, char** argv)
{
	if(argc < 6)
	{
		std::cout << "Usage: <kvalue> <adj list> <lengths file> <estimate file> <max nodes to explore>\n";
		exit(1);
	}
	
	if (string(argv[1]) == "-v") {
		gDebugPrint = true;
		argv++;
		argc--;
	}

	int kmer = atoi(argv[1]);
	std::string adjFile(argv[2]);
	std::string lenFile(argv[3]);
	std::string estFile(argv[4]);
	int maxCost = atoi(argv[5]);
	
	bool preallocVecs = false;
	size_t preallocSize = 0;
	if(argc == 7)
	{
		preallocVecs = true;
		preallocSize = atoi(argv[5]);
	}
	
	std::cout << "Adj file: " << adjFile << " Estimate File: " << estFile << " len file: " << lenFile << " kmer " << kmer << " max cost " << maxCost << " prealloc size " << preallocSize << "\n";
	
	// Create the graph
	SimpleContigGraph* pContigGraph;
	
	if(preallocVecs)
	{
		pContigGraph = new SimpleContigGraph(preallocSize);
	}
	else
	{
		pContigGraph = new SimpleContigGraph();
	}
	
	
	// Load the graph from the adjacency file
	loadGraphFromAdjFile(pContigGraph, lenFile, adjFile);
	
	// try to find paths that match the distance estimates
	generatePathsThroughEstimates(pContigGraph, estFile, kmer, maxCost);
	
	delete pContigGraph;
}

/** Return the allowed error for the given estimate. */
unsigned allowedError(float stddev)
{
	/** The number of standard deviations. */
	const int NUM_SIGMA = 3;

	/**
	 * Additional constant error. The error expected that does not
	 * vary with the number of samples.
	 */
	const unsigned CONSTANT_ERROR = 6;

	return (unsigned)ceilf(NUM_SIGMA * stddev + CONSTANT_ERROR);
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
	std::ofstream outStream("ResolvedPaths.txt");
	SimpleDataCost costFunctor(kmer);

	for (EstimateRecord er; inStream >> er;) {
		for(size_t dirIdx = 0; dirIdx <= 1; ++dirIdx)
		{

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
				
				if(gDebugPrint) std::cout << "Adding Constraint " << iter->nID << " " << nc << "\n";

				constraintMap[iter->nID] = nc;
			}

			if(gDebugPrint)
			{
				std::cout << "Constraints: ";
				SetOperations::printMap(constraintMap);
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
			if (solutions.size() == 0) {
				noPossiblePaths++;
				continue;
			}

			for (SimpleContigGraph::FeasiblePaths::iterator
					solIter = solutions.begin();
					solIter != solutions.end();) {
				if (gDebugPrint) {
					SetOperations::printPath(*solIter);
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
					bool invalid = abs(diff) > buffer;
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
					SetOperations::printPath(*solIter);
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

			if (solutions.size() == 0) {
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
