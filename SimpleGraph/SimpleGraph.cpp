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
void generatePathsThroughEstimates(SimpleContigGraph* pContigGraph, std::string estFileName, int kmer);
void constructContigPath(const SimpleContigGraph::VertexPath& vertexPath, ContigPath& contigPath);
void outputContigPath(std::ofstream& outStream, LinearNumKey refNode, extDirection dir, const ContigPath& contigPath);

//
//
//
int main(int argc, char** argv)
{
	if(argc < 5)
	{
		std::cout << "Usage: <kvalue> <adj list> <lengths file> <estimate file>\n";
		exit(1);
	}
	
	int kmer = atoi(argv[1]);
	std::string adjFile(argv[2]);
	std::string lenFile(argv[3]);
	std::string estFile(argv[4]);
	
	bool preallocVecs = false;
	size_t preallocSize = 0;
	if(argc == 6)
	{
		preallocVecs = true;
		preallocSize = atoi(argv[5]);
	}
	
	std::cout << "Adj file: " << adjFile << " Estimate File: " << estFile << " len file: " << lenFile << " kmer " << kmer << std::endl;
	
	// Create the graph
	SimpleContigGraph* pContigGraph;
	ContigLengthVec contigLens;	
	
	if(preallocVecs)
	{
		std::cout << "preallocating vectors to be " << preallocSize << "\n";
		pContigGraph = new SimpleContigGraph(preallocSize);
		contigLens.reserve(preallocSize);
	}
	else
	{
		pContigGraph = new SimpleContigGraph();
	}
	
	// Load the lengths
	loadContigLengths(lenFile, contigLens);
	
	// Load the graph from the adjacency file
	loadGraphFromAdjFile(pContigGraph, contigLens, adjFile);
	
	// try to find paths that match the distance estimates
	generatePathsThroughEstimates(pContigGraph, estFile, kmer);
	
	delete pContigGraph;
}

//
//
//
void generatePathsThroughEstimates(SimpleContigGraph* pContigGraph, std::string estFileName, int kmer)
{
	int totalAttempted = 0;
	int nopathEnd = 0;
	int uniqueEnd = 0;
	int multiEnd = 0;

	(void)pContigGraph;
	std::ifstream inStream(estFileName.c_str());
	std::ofstream outStream("ResolvedPaths.txt");
	SimpleDataCost costFunctor(kmer);

	while(!inStream.eof() && inStream.peek() != EOF)
	{
		EstimateRecord er;
		readEstimateRecord(inStream, er);
		
		for(size_t dirIdx = 0; dirIdx <= 1; ++dirIdx)
		{

			std::cout << "\n\n****Processing " << er.refID << " dir: " << dirIdx << "****\n\n";
			// generate the reachable set
			SimpleContigGraph::KeyIntMap constraintMap;
			
			const int distanceBuffer = 5;
			for(EstimateVector::iterator iter = er.estimates[dirIdx].begin(); iter != er.estimates[dirIdx].end(); ++iter)
			{
				// Translate the distances produced by the esimator into the coordinate space
				// used by the graph (a translation of m_overlap)
				int translatedDistance = iter->distance + costFunctor.m_overlap;
				constraintMap[iter->nID] = translatedDistance  + distanceBuffer;
			}

			std::cout << "Constraints: ";
			SetOperations::printMap(constraintMap);
			std::cout << "\n";

			// Generate the paths
			SimpleContigGraph::FeasiblePaths solutions;
			
			// The value at which to abort trying to find a path
			// only unique paths are being searched for...
			const int maxNumPaths = 2;
			pContigGraph->findSuperpaths(er.refID, (extDirection)dirIdx, constraintMap, solutions, costFunctor, maxNumPaths);
			
			std::cout << "Solutions: \n";
			for(SimpleContigGraph::FeasiblePaths::iterator solIter = solutions.begin(); solIter != solutions.end(); ++solIter)
			{
				size_t len = pContigGraph->calculatePathLength(*solIter, costFunctor);
				std::cout << len << " ";
				SetOperations::printPath(*solIter);
				std::cout << "\n";
			}
			
			if(solutions.size() == 1)
			{
				printf("Unique solution found!\n");
				bool validPath = true;
				// Calculate the path distance to each node and see if it is within the estimated distance
				SimpleContigGraph::VertexPath sol = solutions.front();
				
				std::map<LinearNumKey, int> distanceMap;
				pContigGraph->makeDistanceMap(sol, costFunctor, distanceMap);
				
				// Filter
				for(EstimateVector::iterator iter = er.estimates[dirIdx].begin(); iter != er.estimates[dirIdx].end(); ++iter)
				{
					std::cout << "Contig " << *iter << std::endl;
					
					// look up in the distance map
					std::map<LinearNumKey, int>::iterator dmIter = distanceMap.find(iter->nID);
					if(dmIter != distanceMap.end())
					{
						// translate distance by -overlap to match coordinate space used by the estimate
						int actualDistance = dmIter->second - costFunctor.m_overlap;
						std::cout << " Actual Dist: " << actualDistance << std::endl;
						int diff = abs(actualDistance - iter->distance);
						std::cout << " diff: " << diff << std::endl;
						
						// Arbitrary cutoff for now
						const int dist_cutoff = 5;
						if(diff > dist_cutoff)
						{
							validPath = false;
							break;
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
					constructContigPath(sol, cPath);
					outputContigPath(outStream, er.refID, (extDirection)dirIdx, cPath);	
				}
				uniqueEnd++;
			}
			else if(solutions.size() > 1)
			{
				multiEnd++;
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
	std::cout << "outting path " << refNode << std::endl;
	outStream << "@ " << refNode << "," << dir << " -> ";
	outStream << contigPath << "\n";
}
