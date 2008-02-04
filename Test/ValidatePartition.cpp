#include <string>
#include "Config.h"
#include "PackedSeqReader.h"
#include "PhaseSpace.h"

int main(int argv, char** argc)
{

	std::string configFile = argc[1];
	
	Config config;
	config.readConfig(configFile);

	int count = 0;	
	for(int x = 0; x < config.getSequenceLength(); x += config.getUnitSize())
		for(int y = 0; y < config.getSequenceLength(); y += config.getUnitSize())
			for(int z = 0; z < config.getSequenceLength(); z += config.getUnitSize())
				for(int w = 0; w < config.getSequenceLength(); w += config.getUnitSize())
				{
					char filename[512];
					sprintf(filename, "%s/x_%d/y_%d/z_%d/w_%d/%s", config.getRootDataDir().c_str(), x, y, z, w, config.getSequenceFilename().c_str());
					
					// open the file
					PackedSeqReader reader(filename);
					

					while(reader.isGood())
					{
						
						PackedSeq seq = reader.ReadSequence();
						printf(">%d\n%s\n", count, seq.decode().c_str());
						/*
						// Check if the sequence is in the correct bucket and the correct size
						if(seq.decode().length() != config.getSequenceLength())
						{
							printf("incorrect sequence length! (%d != %d) for %s\n", seq.decode().length(), config.getSequenceLength(), seq.decode().c_str());
						}
						
						// Check if the sequence is in the correct bin
						Coord4 c = PhaseSpace::SequenceToCoord4(seq);
						
						if(!((c.x >= x && c.x < x + config.getUnitSize()) &&
							(c.y >= y && c.y < y + config.getUnitSize()) &&
							(c.z >= z && c.z < z + config.getUnitSize()) &&
							(c.w >= w && c.w < w + config.getUnitSize())))
						{
							printf("sequence is in incorrect bin coord: (%d %d %d %d) bin start: (%d %d %d %d)\n", c.x, c.y, c.z, c.w, x, y, z, w);
						}
						*/
						count++;
					}
					
					//printf("%d %d %d %d %d\n", x,y,z,w, count);
				}
}
