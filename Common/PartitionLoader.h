#ifndef PARTITIONLOADER_H
#define PARTITIONLOADER_H

#include "CommonDefs.h"
#include "PhaseSpace.h"
#include "PackedSeqReader.h"
#include "Config.h"

// Wrapper class to take in config values and the dimensions of the partition we want
// It finds the files on the system and loads them into the phasespace
class PartitionLoader
{
	public:
		PartitionLoader(const Config* pConfig);
		PhaseSpace* CreateAndLoadPhaseSpace(Coord4 start, Coord4 size);
		static std::string Coord4ToPartitionFile(const Config* pConfig, const Coord4& pos);
		
	private:
		const Config* m_pConfig;
};

#endif
