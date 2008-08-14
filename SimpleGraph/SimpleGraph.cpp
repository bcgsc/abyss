#include <stdio.h>
#include <math.h>
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
#include "VisitAlgorithms.h"
#include "ContigPath.h"

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

bool gDebugPrint = false;

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

//
//
//
void generatePathsThroughEstimates(SimpleContigGraph* pContigGraph, std::string estFileName, int kmer, int maxCost)
{
	int totalAttempted = 0;
	int nopathEnd = 0;
	int uniqueEnd = 0;
	int multiEnd = 0;

	(void)pContigGraph;
	std::ifstream inStream(estFileName.c_str());
	std::ofstream outStream("ResolvedPaths.txt");
	SimpleDataCost costFunctor(kmer);

	// How many standard deviations to look for the estimate
	const int NUM_SIGMA = 3;
	
	while(!inStream.eof() && inStream.peek() != EOF)
	{
		EstimateRecord er;
		readEstimateRecord(inStream, er);
		
		for(size_t dirIdx = 0; dirIdx <= 1; ++dirIdx)
		{

			std::cout << "****Processing " << er.refID << " dir: " << dirIdx << "****\n";
			// generate the reachable set
			SimpleContigGraph::KeyConstraintMap constraintMap;
			
			for(EstimateVector::iterator iter = er.estimates[dirIdx].begin(); iter != er.estimates[dirIdx].end(); ++iter)
			{
				// Translate the distances produced by the esimator into the coordinate space
				// used by the graph (a translation of m_overlap)
				int translatedDistance = iter->distance + costFunctor.m_overlap;
				int distanceBuffer = iter->stdDev * NUM_SIGMA;
				
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
			// only unique paths are being searched for...
			const int maxNumPaths = 200;
			int numVisited = 0;
			
			// The number of nodes to explore before bailing out
			pContigGraph->findSuperpaths(er.refID, (extDirection)dirIdx, constraintMap, solutions, costFunctor, maxNumPaths, maxCost, numVisited);
			
			// Select the best path
			//std::cout << "Computational cost " << numVisited << "\n";
			//std::cout << "Solutions: \n";

			SimpleContigGraph::FeasiblePaths::iterator bestSol = solutions.end();
			int minDiff = 999999;

			if(solutions.size() == 1)
			{
			    for(SimpleContigGraph::FeasiblePaths::iterator solIter = solutions.begin(); solIter != solutions.end(); ++solIter)
			    {
					size_t len = pContigGraph->calculatePathLength(*solIter, costFunctor);
					
					if(gDebugPrint)
					{
						std::cout << len << " ";
						SetOperations::printPath(*solIter);
					}
					
					
					std::map<LinearNumKey, int> distanceMap;
					pContigGraph->makeDistanceMap(*solIter, costFunctor, distanceMap); 
					int sumDiff = 0;
					for(EstimateVector::iterator iter = er.estimates[dirIdx].begin(); iter != er.estimates[dirIdx].end(); ++iter)
					{
					    std::map<LinearNumKey, int>::iterator dmIter = distanceMap.find(iter->nID);
					    if(dmIter != distanceMap.end())
					    {	
							int actualDistance = dmIter->second - costFunctor.m_overlap;
							int diff = actualDistance - iter->distance;
							sumDiff += abs(diff);
					    }
					}
					
					if(sumDiff < minDiff)
					{
						minDiff = sumDiff;
						bestSol = solIter;
					}
					
					if(gDebugPrint) std::cout << " Path diff: " << sumDiff << "\n"; 
			    }
			}
			
			
			
			if(bestSol != solutions.end())
			{
				if(gDebugPrint)
				{
					std::cout << "Selected path with sum diff " << minDiff << std::endl;
					printf("Solution found!\n");
				}
				
				bool validPath = true;
				// Calculate the path distance to each node and see if it is within the estimated distance
				
				std::map<LinearNumKey, int> distanceMap;
				pContigGraph->makeDistanceMap(*bestSol, costFunctor, distanceMap);
				
				// Filter out potentially bad hits
				for(EstimateVector::iterator iter = er.estimates[dirIdx].begin(); iter != er.estimates[dirIdx].end(); ++iter)
				{
					if(gDebugPrint)  std::cout << "Contig " << *iter << "\n";
					
					// look up in the distance map
					std::map<LinearNumKey, int>::iterator dmIter = distanceMap.find(iter->nID);
					if(dmIter != distanceMap.end())
					{
						// translate distance by -overlap to match coordinate space used by the estimate
						int actualDistance = dmIter->second - costFunctor.m_overlap;
						if(gDebugPrint) std::cout << " Actual Dist: " << actualDistance << "\n";
						int diff = actualDistance - iter->distance;
						int buffer = iter->stdDev * NUM_SIGMA;
						if(abs(diff) > buffer)
						{
							validPath = false;
						}
						else
						{
							if(gDebugPrint) std::cout << "diff: " << diff << " n: " << iter->numPairs << " buffer: " << buffer << "\n";
						}
					}
					else
					{
						validPath = false;
					}
				}
				
				// If the path is valid, output it.
				if(validPath)
				{
					ContigPath cPath;
					constructContigPath(*bestSol, cPath);
					outputContigPath(outStream, er.refID, (extDirection)dirIdx, cPath);	
				}
				else
				{
					if(gDebugPrint) printf("Unique path was not valid!\n");
				}
				
				uniqueEnd++;
			}
			else
			{
				nopathEnd++;
			}
			totalAttempted++;
		}
	}
	
	printf("Total paths attempted: %d\n", totalAttempted);
	printf("No paths end: %d\n", nopathEnd);
	printf("Unique end: %d\n", uniqueEnd);
	printf("Multi end: %d\n", multiEnd);
	
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
	std::cout << "Found path for " << refNode << "\n";
	outStream << "@ " << refNode << "," << dir << " -> ";
	outStream << contigPath << "\n";
}
