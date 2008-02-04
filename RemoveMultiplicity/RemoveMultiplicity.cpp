#include <fstream>
#include <set>
#include "RemoveMultiplicity.h"
#include "Reader.h"
#include "PartitionLoader.h"
#include "PackedSeqWriter.h"

int main(int argv, char** argc)
{	
	Reader fileReader;
	
	if(argv < 2)
	{
		printf("usage: Trimmer <config file> <control file>\n");
		exit(1);
	}
	
	// get the command line arguments
	std::string configFile = argc[1];
	std::string controlFile = argc[2];
	
	// Read in the config
	Config config;
	config.readConfig(configFile);
	
	// Read the control file and operate on it	
	ReadControlAndRemove(&config, controlFile);
}

void ReadControlAndRemove(Config* pConfig, std::string controlFile)
{
	
	// open a file, read each line and trim the input bin by the multiplicity
	std::ifstream file(controlFile.c_str());
	
	const int MAX_LINE_LENGTH = 512;
	char buffer[MAX_LINE_LENGTH];

	while(!file.eof() && file.peek() != EOF)
	{
		file.getline(buffer, MAX_LINE_LENGTH);
		Coord4 start;
		
		// iterate over the tokens in the line (space delimited)
		if(sscanf(buffer, "%d %d %d %d", &start.x, &start.y, &start.z, &start.w) != 4)
		{
			printf("invalid control file format\n");
			return;
		}	
		
		// Perform the actual removal
		RemoveMultiplicityByCoord(pConfig, start);
	}
}

void RemoveMultiplicityByCoord(Config* pConfig, Coord4 start)
{

	Coord4 minCoords;
	Coord4 maxCoords;
	
	// get the step size
	int stepSize = pConfig->getUnitSize();
	
	// Compute the coordinate extents
	minCoords.x = start.x;
	maxCoords.x = start.x + stepSize - 1;
	
	minCoords.y = start.y;
	maxCoords.y = start.y + stepSize - 1;
	
	minCoords.z = start.z;
	maxCoords.z = start.z + stepSize - 1;
	
	minCoords.w = start.w;
	maxCoords.w = start.w + stepSize - 1;
	
		
	// Create the partition loader
	PartitionLoader pl(pConfig);
	
	// Create the phase space based on the coordinates
	PhaseSpace* pPS = pl.CreateAndLoadPhaseSpace(minCoords, maxCoords);
	
	// Create the output file
	std::string outFile = PartitionLoader::Coord4ToPartitionFile(pConfig, start);
	std::string tempFile = outFile + ".temp";
	printf("Opening %s for write\n", tempFile.c_str());
	
	// Create the writer object
	PackedSeqWriter writer(tempFile.c_str(), pConfig->getSequenceLength());
	
	// Output the sequences
	int count = 0;
	int passed = 0;
	
	// output all the unique sequences in the phase space
	for(int x = 0; x < stepSize; x++)
		for(int y = 0; y < stepSize; y++)
			for(int z = 0; z < stepSize; z++)
				for(int w = 0; w < stepSize; w++)
				{
					Coord4 c;
					c.x = start.x + x;
					c.y = start.y + y;
					c.z = start.z + z;
					c.w = start.w + w;
					
					// The set of sequences already outputted
					std::set<PackedSeq> outputtedSeqs;
					
					// Iterators
					PhaseSpaceBinIter startIter = pPS->getStartIter(c);
					PhaseSpaceBinIter endIter = pPS->getEndIter(c);
					
					// Loop over all the sequences in the bin
					for(PhaseSpaceBinIter iter = startIter; iter != endIter; iter++)
					{
						const PackedSeq& seq = *iter;
						if(outputtedSeqs.count(*iter) == 0 && outputtedSeqs.count(reverseComplement(*iter)) == 0)
						{					
							writer.WriteSequence(*iter);
							outputtedSeqs.insert(*iter);
							passed++;
						}
						//printf("%s %d\n", iter->first.decode().c_str(), iter->second);
						count++;
					}	
					//return;				
				}
	
	printf("%d/%d passed multiplicity\n", passed, count);
	
	// delete the phase space
	delete pPS;
	pPS = 0;
}
