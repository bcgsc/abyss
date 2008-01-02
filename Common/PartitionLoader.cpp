#include <set>
#include "PartitionLoader.h"

PartitionLoader::PartitionLoader(const Config* pConfig) : m_pConfig(pConfig)
{
	
}


PhaseSpace* PartitionLoader::CreateAndLoadPhaseSpace(Coord4 start, Coord4 size)
{
	// Create the phase space
	PhaseSpace* pPS = new PhaseSpace(m_pConfig->getSequenceLength(), start, size);

	// Calculate the files that need to be loaded to determine this slice of 4D space

	// the filenames will be added to this set which guarentees each file will only be read once
	std::set<std::string> fileSet;

	for(int x = start.x; x < start.x + size.x; x++)
		for(int y = start.y; y < start.y + size.y; y++)
			for(int z = start.z; z < start.z + size.z; z++)
				for(int w = start.w; w < start.w + size.w; w++)
				{
					Coord4 pos = {x,y,z,w};
					std::string file = Coord4ToPartitionFile(pos);
					fileSet.insert(file);
				}

	// Read in the files needed to load up the specified hypercube slice
	for(std::set<std::string>::const_iterator iter = fileSet.begin(); iter != fileSet.end(); iter++)
	{
		printf("file: %s\n", iter->c_str());
		PackedSeqReader reader(iter->c_str());
		while(reader.isGood())
		{
			// this call allocates memory
			PackedSeq* s = reader.ReadSequence();
			
			// this call makes a copy of the packed seq
			pPS->addSequence(*s);
			
			// delete the first one
			delete s;
			s = 0;
		}
	}
	
	return pPS;	
}

// Function that computes the correct sequence file to load from a coord4 position
std::string PartitionLoader::Coord4ToPartitionFile(const Coord4& pos) const
{
	int stride = m_pConfig->getUnitSize();
	int x = (pos.x / stride) * stride;
	int y = (pos.y / stride) * stride;
	int z = (pos.z / stride) * stride;
	int w = (pos.w / stride) * stride;	

	printf("(%d %d %d %d) maps to (%d %d %d %d)\n", pos.x, pos.y, pos.z, pos.w, x,y,z,w);
	
	char buffer[512];
	sprintf(buffer, "%s/x_%d/y_%d/z_%d/w_%d/%s", m_pConfig->getRootDataDir().c_str(), x,y,z,w,m_pConfig->getSequenceFilename().c_str());
	std::string s(buffer);
	return s;
}
