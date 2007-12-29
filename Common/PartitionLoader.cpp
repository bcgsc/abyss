#include "PartitionLoader.h"

PartitionLoader::PartitionLoader(const Config* pConfig) : m_pConfig(pConfig)
{
	
}


PhaseSpace* PartitionLoader::CreateAndLoadPhaseSpace(Coord4 start, Coord4 size)
{
	// Create the phase space
	PhaseSpace* pPS = new PhaseSpace(m_pConfig->getSequenceLength(), start, size);

	// Calculate the files that need to be loaded to determine this slice of 4D space
	

	// Read in the files needed to load up the specified hypercube slice

	return pPS;	
}

// Function that computes the correct sequence file to load from a coord4 position
std::string PartitionLoader::Coord4ToPartitionFile(const Coord4& pos) const
{
	int stride = m_pConfig->getPartitionStep();
	int x = (pos.x / stride) * stride; // this looks it should do nothing but it will truncate the value to the lowest multiple of stride...the bin we want
	int y = (pos.y / stride) * stride;
	int z = (pos.z / stride) * stride;
	int w = (pos.w / stride) * stride;	

	printf("(%d %d %d %d) maps to (%d %d %d %d)\n", x,y,z,w);
}
