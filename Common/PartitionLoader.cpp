#include <set>
#include "PartitionLoader.h"

PartitionLoader::PartitionLoader(const Config* pConfig) : m_pConfig(pConfig)
{
	
}


PhaseSpace* PartitionLoader::CreateAndLoadPhaseSpace(Coord4 center, Coord4 size)
{


	// Calculate the coordinate extents
	Coord4 minCoords;
	Coord4 maxCoords;
	minCoords.x = max(center.x - size.x, 0);
	maxCoords.x = min(center.x + size.x + m_pConfig->getUnitSize(), m_pConfig->getSequenceLength());
	
	minCoords.y = max(center.y - size.y, 0);
	maxCoords.y = min(center.y + size.y + m_pConfig->getUnitSize(), m_pConfig->getSequenceLength());
	
	minCoords.z = max(center.z - size.z, 0);
	maxCoords.z = min(center.z + size.z + m_pConfig->getUnitSize(), m_pConfig->getSequenceLength());
	
	minCoords.w = max(center.w - size.w, 0);
	maxCoords.w = min(center.w + size.w + m_pConfig->getUnitSize(), m_pConfig->getSequenceLength());
	
	printf("mins: (%d, %d, %d, %d)\n", minCoords.x, minCoords.y, minCoords.z, minCoords.w);
	printf("maxs: (%d, %d, %d, %d)\n", maxCoords.x, maxCoords.y, maxCoords.z, maxCoords.w);
	
	// Create the phase space
	PhaseSpace* pPS = new PhaseSpace(m_pConfig->getSequenceLength(), minCoords, maxCoords);
		
	// Calculate the files that need to be loaded to determine this slice of 4D space
	std::set<std::string> fileSet;	
	for(int x = minCoords.x; x < maxCoords.x; x++)
		for(int y = minCoords.y; y < maxCoords.y; y++)
			for(int z = minCoords.z; z < maxCoords.z; z++)
				for(int w = minCoords.w; w < maxCoords.w; w++)
				{
					Coord4 pos = {x,y,z,w};
					std::string file = Coord4ToPartitionFile(m_pConfig, pos);
					fileSet.insert(file);
					//printf("coord (%d, %d, %d, %d)\n", x,y,z,w);
				}
				
	int count = 0;
	// Read in the files needed to load up the specified hypercube slice
	for(std::set<std::string>::const_iterator iter = fileSet.begin(); iter != fileSet.end(); iter++)
	{
		printf("file: %s\n", iter->c_str());
		PackedSeqReader reader(iter->c_str());
		while(reader.isGood())
		{
			// this call allocates memory
			PackedSeq s = reader.ReadSequence();
			
			// this call makes a copy of the packed seq
			pPS->addSequence(s);
			count++;
		}
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

int PartitionLoader::min(const int& n1, const int& n2) const
{
	return (n1 < n2) ? n1 : n2;	
}

int PartitionLoader::max(const int& n1, const int& n2) const
{
	return (n1 > n2) ? n1 : n2;	
}
