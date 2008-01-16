#include <set>
#include "PartitionLoader.h"

PartitionLoader::PartitionLoader(const Config* pConfig) : m_pConfig(pConfig)
{
	
}

// Coordinates are inclusive
PhaseSpace* PartitionLoader::CreateAndLoadPhaseSpace(Coord4 minCoord, Coord4 maxCoord)
{

	assert(minCoord.x < maxCoord.x && minCoord.y < maxCoord.y && minCoord.z < maxCoord.z && minCoord.w < maxCoord.w);
	
	// Compute the starting bin coordinates. They are the lowest coordinates that are a multiple of the unitsize
	Coord4 start;
	start.x = (minCoord.x / m_pConfig->getUnitSize()) * m_pConfig->getUnitSize();
	start.y = (minCoord.y / m_pConfig->getUnitSize()) * m_pConfig->getUnitSize();
	start.z = (minCoord.z / m_pConfig->getUnitSize()) * m_pConfig->getUnitSize();
	start.w = (minCoord.w / m_pConfig->getUnitSize()) * m_pConfig->getUnitSize();
	
	// Compute the end bin coordinates. 
	Coord4 end;
	end.x = (maxCoord.x / m_pConfig->getUnitSize()) * m_pConfig->getUnitSize();
	end.y = (maxCoord.y / m_pConfig->getUnitSize()) * m_pConfig->getUnitSize();
	end.z = (maxCoord.z / m_pConfig->getUnitSize()) * m_pConfig->getUnitSize();
	end.w = (maxCoord.w / m_pConfig->getUnitSize()) * m_pConfig->getUnitSize();
		
	printf("start load bin: (%d, %d, %d, %d)\n", start.x, start.y, start.z, start.w);
	printf("end load bin: (%d, %d, %d, %d)\n", end.x, end.y, end.z, end.w);
	
	// Create the phase space
	PhaseSpace* pPS = new PhaseSpace(m_pConfig->getSequenceLength(), minCoord, maxCoord);
	
	int count = 0;
	
	// Load up all the files required
	for(int x = start.x; x <= end.x; x += m_pConfig->getUnitSize())
		for(int y = start.y; y <= end.y; y += m_pConfig->getUnitSize())
			for(int z = start.z; z <= end.z; z += m_pConfig->getUnitSize())
				for(int w = start.w; w <= end.w; w += m_pConfig->getUnitSize())
				{
					Coord4 pos = {x,y,z,w};
					std::string file = Coord4ToPartitionFile(m_pConfig, pos);
					printf("loading...%s\n", file.c_str());
					
					// read in all the sequences
					PackedSeqReader reader(file.c_str());
					while(reader.isGood())
					{
						PackedSeq s = reader.ReadSequence();
						//printf("seq: %s\n", s.decode().c_str());
						pPS->addSequence(s);
						count++;
					}
					
					// finalize the bins
					
					// make sure the coordinates do not go out of bounds (we could be just reading in a subset of the coordinates
					Coord4 startBin;
					startBin.x = max(pos.x, minCoord.x);
					startBin.y = max(pos.y, minCoord.y);
					startBin.z = max(pos.z, minCoord.z);
					startBin.w = max(pos.w, minCoord.w);					
					
					Coord4 endBin;
					endBin.x = min(pos.x + m_pConfig->getUnitSize() - 1, maxCoord.x);
					endBin.y = min(pos.y + m_pConfig->getUnitSize() - 1, maxCoord.y);
					endBin.z = min(pos.z + m_pConfig->getUnitSize() - 1, maxCoord.z);
					endBin.w = min(pos.w + m_pConfig->getUnitSize() - 1, maxCoord.w);
					pPS->finalizeBins(startBin, endBin);
				}
	
	printf("loaded %d sequences\n", count); 
	return pPS;	
}

// Function that computes the correct sequence file to load from a coord4 position
std::string PartitionLoader::Coord4ToPartitionFile(const Config* pConfig, const Coord4& pos)
{
	int stride = pConfig->getUnitSize();
	int x = (pos.x / stride) * stride;
	int y = (pos.y / stride) * stride;
	int z = (pos.z / stride) * stride;
	int w = (pos.w / stride) * stride;	

	//printf("(%d %d %d %d) maps to (%d %d %d %d)\n", pos.x, pos.y, pos.z, pos.w, x,y,z,w);
	
	char buffer[512];
	sprintf(buffer, "%s/x_%d/y_%d/z_%d/w_%d/%s", pConfig->getRootDataDir().c_str(), x,y,z,w,pConfig->getSequenceFilename().c_str());
	std::string s(buffer);
	return s;
}

