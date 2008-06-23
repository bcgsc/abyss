#include <stdio.h>
#include <math.h>
#include <iostream>
#include "SequenceCollectionHash.h"
#include "AssemblyAlgorithms.h"
#include "Options.h"
#include "FastaReader.h"
#include "Stats.h"
#include "AlignExtractor.h"


typedef std::map<ContigID, int> ContigLengthMap;

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
void loadContigLengths(std::string contigLenFile, ContigLengthMap& lengthMap);
int lookupLength(const ContigLengthMap& lengthMap, const ContigID& id);
void processContigs(std::string alignFile, const ContigLengthMap& lengthMap, const PDF& pdf);
PDF loadPDF(std::string distCountFile, const int limit);
void writeRecord(std::ofstream& stream, ContigID id, int distance, int n);

// Go through a list of pairings and provide a maximum likelihood estimate of the distance
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
		ContigID refContigID = currPairs.front().refRec.contigID;
		
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
					PairedData& pd = dataMap[iter->pairRec.contigID];
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
				
				// Estimate the distance
				int distance = estimateDistance(refContigLength, pairContigLength, dirIdx, pdIter->second, pdf);
				
				// write the record to file
				writeRecord(outFile, pairID, distance, pdIter->second.pairVec.size());
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
			const int read_len = 36;
			apIter->pairRec.position = (pairLen - (apIter->pairRec.position + read_len)); 
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
			int refTransPos = apIter->refRec.position + refOffset;
			int pairTransPos = apIter->pairRec.position + pairOffset;
			
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
	return dist;
}

void writeRecord(std::ofstream& stream, ContigID id, int distance, int n)
{
	stream << id << "," << distance << "," << n << " ";
}

void loadContigLengths(std::string contigLenFile, ContigLengthMap& lengthMap)
{
	ifstream contigLenStream(contigLenFile.c_str());
	while(!contigLenStream.eof() && contigLenStream.peek() != EOF)
	{
		ContigID id;
		int len;
		contigLenStream >> id >> len;
		lengthMap[id] = len;

	}
	contigLenStream.close();
}

int lookupLength(const ContigLengthMap& lengthMap, const ContigID& id)
{
	ContigLengthMap::const_iterator iter = lengthMap.find(id);
	assert(iter != lengthMap.end());
	return iter->second;
}

PDF loadPDF(std::string distCountFile, const int limit)
{
	Histogram hist;
	ifstream distFile(distCountFile.c_str());
	while(!distFile.eof() && distFile.peek() != EOF)
	{
		int value;
		int count;
		distFile >> value;
		distFile >> count;

		if(value < limit)
		{
			//std::cout << "adding " << value << " : " << count << std::endl;
			hist.addMultiplePoints(value, count);
		}
	} 

	PDF pdf(hist);
	return pdf;
}




