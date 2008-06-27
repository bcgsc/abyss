#include <stdio.h>
#include <math.h>
#include <iostream>
#include "SequenceCollectionHash.h"
#include "AssemblyAlgorithms.h"
#include "Options.h"
#include "FastaReader.h"
#include "Stats.h"
#include "AlignExtractor.h"
#include "PairUtils.h"

typedef std::vector<int> IntVector;
struct PairedData
{
	PairedData()
	{
		compCount[0] = 0;
		compCount[1] = 0;
	}

	AlignPairVec pairVec;
	int compCount[2];
};

typedef std::map<ContigID, PairedData> PairDataMap;

// FUNCTIONS
int estimateDistance(int refLen, int pairLen, size_t dirIdx, PairedData& pairData, const PDF& pdf);
void processContigs(std::string alignFile, const ContigLengthMap& lengthMap, const PDF& pdf);

/*
// Go through a list of pairings and provide a maximum likelihood estimate of the distance
int main(int argc, char** argv)
{
	(void)argc;
	std::string distanceCountFile(argv[1]);
	std::string distanceListFile(argv[2]);
	PDF empiricalPDF = loadPDF(distanceCountFile, 350);
	
	std::ifstream inFile(distanceListFile.c_str());
	const int num_samples = atoi(argv[3]);
	
	IntVector distances;
	
	for(int i = 0; i < num_samples; ++i)
	{
		int d;
		inFile >> d;
		if(d < 200)
			distances.push_back(d);
	}
	
	KSTestCont(distances, empiricalPDF);
	
	return 1;
}
*/

int main(int argc, char** argv)
{
	if(argc < 4)
	{
		std::cout << "Usage: <alignFile> <length file> <distance count file>\n";
		exit(1);
	}
	
	std::string alignFile(argv[1]);
	std::string contigLengthFile(argv[2]);
	std::string distanceCountFile(argv[3]);

	std::cout << "Align File: " << alignFile << "len file: " << contigLengthFile << " distance file: " << distanceCountFile << std::endl;
	std::cout << "WARNING UNHARDCODE READ LENGTH\n";
	PDF empiricalPDF = loadPDF(distanceCountFile, 350);

	
	// Load the length map
	ContigLengthMap contigLens;	
	loadContigLengths(contigLengthFile, contigLens);

	// Estimate the distances between contigs, one at a time
	processContigs(alignFile, contigLens, empiricalPDF);
	
	return 1;
} 

void processContigs(std::string alignFile, const ContigLengthMap& lengthMap, const PDF& pdf)
{
	(void)pdf;
	AlignExtractor extractor(alignFile);
	
	// open the output file
	std::ofstream outFile("EstimatedLinks.txt");
	
	int count = 0;
	//Extract the align records from the file, one contig's worth at a time
	bool stop = false;
	
	while(!stop)
	{
		
		AlignPairVec currPairs;
		stop = extractor.extractContigAlignments(currPairs);

		assert(currPairs.size() > 0);
		ContigID refContigID = currPairs.front().refRec.contig;
		
		std::cout << "Ref ctg " << refContigID << "\n";
		// Only process contigs that are a reasonable length
		int refLength = lookupLength(lengthMap, refContigID);
		if(refLength < 100)
		{
			continue;
		}
		
		// Write the first field to the file
		outFile << refContigID << " : ";
		
		//std::cout << "Contig " << refContigID << " has " << currPairs.size() << " alignments\n";

		// Seperate the pairings by direction (pairs aligning in the same comp as the contig
		// are sense pairs) and by the contig they align to
		for(size_t dirIdx = 0; dirIdx <= 1; ++dirIdx)
		{
			// If this is the second direction, write a seperator
			if(dirIdx == 1)
			{
				outFile << " | ";
			}
			
			PairDataMap dataMap;
			for(AlignPairVec::iterator iter = currPairs.begin(); iter != currPairs.end(); ++iter)
			{
				if(iter->refRec.isRC == (bool)dirIdx)
				{
					PairedData& pd = dataMap[iter->pairRec.contig];
					pd.pairVec.push_back(*iter);
					size_t compIdx = (size_t)iter->pairRec.isRC;
					assert(compIdx < 2);
					pd.compCount[compIdx]++;
				}
			}

			// For each contig that is paired, compute the distance
			for(PairDataMap::iterator pdIter = dataMap.begin(); pdIter != dataMap.end(); ++pdIter)
			{
				ContigID pairID = pdIter->first;
				// get the contig lengths
				int refContigLength = lookupLength(lengthMap, refContigID);
				int pairContigLength = lookupLength(lengthMap, pairID);
				
				// Check if the pairs are in a valid orientation
				if(pdIter->second.compCount[0] > 0 && pdIter->second.compCount[1] > 0)
				{
					// The pairs are inconsistent, do not use this contig
					continue;
				}
				
				const size_t number_of_pairs_threshold = 30;
				
				if(pdIter->second.pairVec.size() > number_of_pairs_threshold)
				{
					// Estimate the distance
					int distance = estimateDistance(refContigLength, pairContigLength, dirIdx, pdIter->second, pdf);
					
					Estimate est;
					est.cID = pairID;
					est.distance = distance;
					est.numPairs = pdIter->second.pairVec.size();
					
					// write the record to file
					outFile << est << " "; 
				}
				//std::cout << "Est dist: " << dist << " ratio " << ratio << std::endl;
			}
		}
		outFile << "\n";
		count++;
		if(count % 10000 == 0)
		{
			std::cout << "Processed " << count << " contigs\n";
		}
	}
	
	outFile.close();
}

// Estimate the distances between the contigs
int estimateDistance(int refLen, int pairLen, size_t dirIdx, PairedData& pairData, const PDF& pdf)
{
	// Determine the relative orientation of the contigs
	// As pairs are orientated in opposite (reverse comp) direction, the alignments are in the same
	// orientation if the pairs aligned in the opposite orientation
	bool sameOrientation = false;
	if(pairData.compCount[1 - dirIdx] > 0)
	{
		// the contigs have the correct orientation
		sameOrientation = true;
	}
	else
	{
		sameOrientation = false;
	}

	// Calculate the distance list for this contig
	// The provisional distances are calculated as if the contigs overlap perfectly by k-1 bases
	// The maximum likelihood estimate will refine this

	// Setup the offsets
	if(!sameOrientation)
	{
		// Flip all the positions of the pair aligns
		for(AlignPairVec::iterator apIter = pairData.pairVec.begin(); apIter != pairData.pairVec.end(); ++apIter)
		{
			apIter->pairRec.start = (pairLen - (apIter->pairRec.start + apIter->pairRec.length)); 
		}
		
	}

	int refOffset = 0;
	int pairOffset = 0;
	
	if(dirIdx == 0)
	{
		// refContig is on the left, offset pairContig by the length of refContig
		pairOffset = refLen;
	}
	else
	{
		// pairContig is on the left, offset refContig by the length of pairContig
		refOffset = pairLen;
	}	
	
	IntVector distanceList;
	for(AlignPairVec::iterator apIter = pairData.pairVec.begin(); apIter != pairData.pairVec.end(); ++apIter)
	{
			int distance;
			int refTransPos = apIter->refRec.start + refOffset;
			int pairTransPos = apIter->pairRec.start + pairOffset;
			
			if(refTransPos < pairTransPos)
			{
				distance = 	pairTransPos - refTransPos;
			}
			else
			{
				distance = 	refTransPos - pairTransPos;
			}
			
			distanceList.push_back(distance);
			//std::cout << "Distance: " << distance << std::endl;
	}
	// Perform the max-likelihood est
	double ratio;
	int dist = maxLikelihoodEst(-25, 200, distanceList, pdf, ratio);
	
	/*
	printf("Num Pairs: %d\n", distanceList.size());
	// Perform a KS-test
	IntVector transDistances;
	for(IntVector::iterator transIter = distanceList.begin(); transIter != distanceList.end(); ++transIter)
	{
		transDistances.push_back(*transIter + dist);
	}
	KSTestCont(transDistances, pdf);
	*/
	
	
	return dist;
}

