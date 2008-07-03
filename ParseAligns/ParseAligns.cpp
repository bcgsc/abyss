#include <stdio.h>

#include <iostream>
#include <fstream>
#include <string>
#include "Aligner.h"
#include "PairUtils.h"

// TYPEDEFS
typedef std::map<std::string, AlignmentVector> ReadAlignMap;


// FUNCTIONS
std::string makePairID(std::string refID);
void readAlignmentsFile(std::string filename, ReadAlignMap& alignTable);

int main(int argc, char* const* argv)
{

	if(argc < 2)
	{
		std::cout << "Usage: <list of files to parse>\n";
		exit(1);
	}
	
	int numFiles = argc - 1;
	
	std::ofstream pairedAlignFile("PairAligns.txt");
	std::ofstream distanceList("DistanceList.txt");	
	
	ReadAlignMap alignTable;
	
	for(int fileIdx = 0; fileIdx < numFiles; ++fileIdx)
	{
		std::string file(argv[fileIdx+1]);
		readAlignmentsFile(file, alignTable);
	}
				
	int numDifferent = 0;
	int numSame = 0;
	int numMissed = 0;
	int numMulti = 0;
	
	// parse the alignment table
	for(ReadAlignMap::iterator iter = alignTable.begin(); iter != alignTable.end(); ++iter)
	{
		std::string currID = iter->first;
		std::string pairID = makePairID(currID);
		
		// Find the pair align
		ReadAlignMap::iterator pairIter = alignTable.find(pairID);
		if(pairIter != alignTable.end())
		{
			if(iter->second.size() == 1 && pairIter->second.size() == 1)
			{
				// Alignments are unique and usable
				Alignment refAlign = iter->second.front();
				Alignment pairAlign = pairIter->second.front();
				
				// Are they on the same contig?
				if(refAlign.contig == pairAlign.contig)
				{
					// Calculate the distance between the reads
					int dist = abs(refAlign.start -  pairAlign.start);
					
					// Print it to the file
					distanceList << dist << "\n";
					numSame++;
				}
				else
				{
					// Print the alignment and the swapped alignment
					pairedAlignFile << currID << " " << refAlign << " " << pairID << " " << pairAlign << "\n";
					pairedAlignFile << pairID << " " << pairAlign << " " << currID << " " << refAlign << "\n";
					numDifferent++;
				}					
			}
			else
			{
				numMulti++;
			}

			
			// Erase the pair as its not needed (faster to mark it as invalid?)
			alignTable.erase(pairIter);
		}
		else
		{
			numMissed++;
		}
	}
	
	printf("Num unmatched: %d Num same contig: %d Num diff contig: %d Num multi: %d\n", numMissed, numSame, numDifferent, numMulti);
	
	pairedAlignFile.close();
	distanceList.close();
}

// Read in the alignments file into the table
void readAlignmentsFile(std::string filename, ReadAlignMap& alignTable)
{
	(void)alignTable;
	std::cout << "Reading " << filename << std::endl;
	std::ifstream fileStream(filename.c_str());
	std::string line;
	while(getline(fileStream,line))
	{
		// Parse the id
		std::string readID;
		
		std::stringstream ss(line);
		
		ss >> readID;
		
		// parse the alignments
		Alignment ali;
		while(ss >> ali)
		{
			alignTable[readID].push_back(ali);
		}
	}
	
	fileStream.close();
}

// Change the id into the id of its pair
std::string makePairID(std::string refID)
{
	std::string pairID = refID;
	// Change the last character
	size_t lastIdx = pairID.size() - 1;
	pairID[lastIdx] = (pairID[lastIdx] == '1') ? '2' : '1';
	return pairID;
}
