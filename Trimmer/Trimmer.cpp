#include <fstream>
#include <set>
#include "Trimmer.h"
#include "PackedSeqWriter.h"

int main(int argv, char** argc)
{	
	if(argv < 3)
	{
		printf("usage: Trimmer <config file> <control file>\n");
		exit(1);
	}
	
	std::string configFile = argc[1];
	std::string controlFile = argc[2];
	
	// Read in the config
	Config config;
	config.readConfig(configFile);
	
	// Read the control file and operate on it	
	ReadControlAndTrim(&config, controlFile);
}

void ReadControlAndTrim(Config* pConfig, std::string controlFile)
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
		TrimByCoordinate(pConfig, start);
	}
}

void TrimByCoordinate(Config* pConfig, Coord4 trimCoord)
{
	// Calculate the coordinate extents, all the sequences in this range will be loaded. these coordinates are INCLUSIVE
	Coord4 minCoords;
	Coord4 maxCoords;
	
	// get the step size
	int stepSize = pConfig->getUnitSize();
	
	// Compute the coordinate extents
	minCoords.x = max(trimCoord.x - 2, 0);
	maxCoords.x = min(trimCoord.x + stepSize + 2, pConfig->getSequenceLength() -1);
	
	minCoords.y = max(trimCoord.y - 1, 0);
	maxCoords.y = min(trimCoord.y + stepSize, pConfig->getSequenceLength() - 1);
	
	minCoords.z = max(trimCoord.z - 1, 0);
	maxCoords.z = min(trimCoord.z + stepSize, pConfig->getSequenceLength() - 1);
	
	minCoords.w = max(trimCoord.w - 1, 0);
	maxCoords.w = min(trimCoord.w + stepSize, pConfig->getSequenceLength() - 1);	

	printf("trim coord: (%d, %d, %d, %d)\n", trimCoord.x, trimCoord.y, trimCoord.z, trimCoord.w);
	printf("min coord: (%d, %d, %d, %d)\n", minCoords.x, minCoords.y, minCoords.z, minCoords.w);
	printf("max coord: (%d, %d, %d, %d)\n", maxCoords.x, maxCoords.y, maxCoords.z, maxCoords.w);
	
	PartitionLoader pl(pConfig);
	PhaseSpace* pPS = pl.CreateAndLoadPhaseSpace(minCoords, maxCoords);	
	
	// Perform the trimming
	// output all the unique sequences in the phase space
	int passedCheck = 0;
	int passedTrimming = 0;
	int count = 0;
	
	// Create the output file
	std::string outFile = PartitionLoader::Coord4ToPartitionFile(pConfig, trimCoord);
	std::string tempFile = outFile + pConfig->getTempFileExtension();	
	
	PackedSeqWriter writer(tempFile.c_str(), pConfig->getSequenceLength());
	printf("opened %s for write\n", tempFile.c_str());
	
	int maxX = min(stepSize, pConfig->getSequenceLength() - trimCoord.x);
	int maxY = min(stepSize, pConfig->getSequenceLength() - trimCoord.y);
	int maxZ = min(stepSize, pConfig->getSequenceLength() - trimCoord.z);
	int maxW = min(stepSize, pConfig->getSequenceLength() - trimCoord.w);
	 
	for(int x = 0; x < maxX; x++)
		for(int y = 0; y < maxY; y++)
			for(int z = 0; z < maxZ; z++)
				for(int w = 0; w < maxW; w++)
				{
					Coord4 c;
					c.x = trimCoord.x + x;
					c.y = trimCoord.y + y;
					c.z = trimCoord.z + z;
					c.w = trimCoord.w + w;
					printf("trimmingcoord [%d %d %d %d]\n", c.x, c.y, c.z, c.w);
					
					for(PhaseSpaceBinIter iter = pPS->getStartIter(c); iter != pPS->getEndIter(c); iter++)
					{
						//printf("%d %s\n", count, iter->decode().c_str());
						//if(pPS->hasChild(*iter) && pPS->hasParent(*iter))
						if(pPS->checkForSequence(*iter))
						{
							passedCheck++;
						}
						
						//printf("trimming...%s\n", iter->decode().c_str());
						
						if(pPS->hasChild(*iter) && pPS->hasParent(*iter))
						{
							//printf("	->passed trimming\n");
							writer.WriteSequence(*iter);
							passedTrimming++;
						}
						count++;
					}					
				}
				
	printf("%d of %d passed check\n", passedCheck, count);
	printf("%d of %d passed trimming\n", passedTrimming, count);
	
	delete pPS;
	pPS = 0;
}
