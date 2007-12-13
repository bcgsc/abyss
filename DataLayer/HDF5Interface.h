#ifndef HDF5INTERFACE_H
#define HDF5FINTERFACE_H

#include <stdio.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#include "PhaseSpace.h"

#define SEQLEN 37
#define NFIELDS  (hsize_t)  2

struct SeqRecord
{
	int id;
	char seq[SEQLEN];
};


#define ERR(e) {printf("NetCDF error: %s\n", nc_strerror(e));}

typedef std::vector<int> Count1D;
typedef std::vector<Count1D> Count2D;
typedef std::vector<Count2D> Count3D;
typedef std::vector<Count3D> Count4D;


class HDF5Interface
{
	public:
		// Constructor opens the file
		HDF5Interface(const char* filename, FileMode mode);
		
		// create a new file and set up all the variables needed to write to it
		void openNewFileForWrite(const char* filename);
		void openFileForRead(const char* filename);
		
		// read a group of sequences
		void readSequencesAtCoordinate(Coord4 pos, Coord4 size);
		
		// write a group of sequences
		void writeSequence(Sequence& seq);
		
		// Destructor closes the file
		~HDF5Interface();


	private:
		static const int NDIMS = 6;
		int dimids[NDIMS];
		int m_seqVarID;
		int m_fileID;
		FileMode m_mode;
		Count4D m_recordCounts;
		
		// table fields
		size_t m_DSTSize;
		size_t m_DSTOffset[NFIELDS];
		size_t m_DSTSizes[NFIELDS];
	
};

#endif
