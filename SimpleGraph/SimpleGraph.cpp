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

// Typedefs
typedef std::vector<NumericID> NumericIDVector;
typedef std::vector<NumericIDVector > NumericIDVectorVector;

// FUNCTORS
struct SimpleDataCost
{
	SimpleDataCost(size_t kmer) : m_overlap(kmer - 1) { }
	
	size_t cost(const SimpleContigData& data) { return data.length - m_overlap; } 
	size_t m_overlap;
};


// FUNCTIONS
void generatePathsThroughEstimates(SimpleContigGraph* pContigGraph, std::string estFileName, int kmer);
void outputPath(std::ofstream& outStream, NumericID refNode, extDirection dir, const SimpleContigGraph::VertexPath& path);

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
	std::cout << "Adj file: " << adjFile << " Estimate File: " << estFile << " len file: " << lenFile << " kmer " << kmer << std::endl;
	
	// Create the graph
	SimpleContigGraph* pContigGraph = new SimpleContigGraph;
	
	// Load the lengths
	ContigLengthMap contigLens;	
	loadContigLengths(lenFile, contigLens);
	
	// Load the graph from the adjacency file
	loadGraphFromAdjFile(pContigGraph, contigLens, adjFile);
	
	// try to find paths that match the distance estimates
	generatePathsThroughEstimates(pContigGraph, estFile, kmer);
}

//
//
//
void generatePathsThroughEstimates(SimpleContigGraph* pContigGraph, std::string estFileName, int kmer)
{
	(void)pContigGraph;
	std::ifstream inStream(estFileName.c_str());
	std::ofstream outStream("ResolvedPaths.txt");

	while(!inStream.eof() && inStream.peek() != EOF)
	{
		EstimateRecord er;
		readEstimateRecord(inStream, er);
		
		for(size_t dirIdx = 0; dirIdx <= 1; ++dirIdx)
		{
			std::cout << "\n\n****Processing " << er.refID << " dir: " << dirIdx << "****\n\n";
			// generate the reachable set
			std::set<NumericID> idSet;
			
			std::cout << "Reachable: { ";
			for(EstimateVector::iterator iter = er.estimates[dirIdx].begin(); iter != er.estimates[dirIdx].end(); ++iter)
			{
				std::cout << *iter << " ";
				idSet.insert(iter->nID);
			}
			
			std::cout << "}\n";
			
			// Generate the paths
			SimpleDataCost costFunctor(kmer);
			SimpleContigGraph::FeasiblePaths solutions;
			
			// The value at which to abort trying to find a path
			// only unique paths are being searched for...
			const int maxNumPaths = 2;
			pContigGraph->findSuperpaths(er.refID, (extDirection)dirIdx, idSet, solutions, costFunctor, maxNumPaths);
			
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
				
				std::map<NumericID, int> distanceMap;
				pContigGraph->makeDistanceMap(sol, costFunctor, distanceMap);
				
				// Filter
				for(EstimateVector::iterator iter = er.estimates[dirIdx].begin(); iter != er.estimates[dirIdx].end(); ++iter)
				{
					std::cout << "Contig " << *iter << std::endl;
					
					// look up in the distance map
					std::map<NumericID, int>::iterator dmIter = distanceMap.find(iter->nID);
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
					outputPath(outStream, er.refID, (extDirection)dirIdx, sol);
				}				
				

			}
			
		}
	}
	
	inStream.close();
	outStream.close();
}

void outputPath(std::ofstream& outStream, NumericID refNode, extDirection dir, const SimpleContigGraph::VertexPath& path)
{
	ContigID realID = convertNumericIDToContigID(refNode);
	outStream << "@ " << realID << "," << dir << " -> ";
	
	for(SimpleContigGraph::VertexPath::const_iterator pIter = path.begin(); pIter != path.end(); ++pIter)
	{
		ContigID id = convertNumericIDToContigID(pIter->key);
		MergeNode mn = {id, pIter->isRC};
		outStream << mn << " ";
	}

	outStream << "\n";
}
