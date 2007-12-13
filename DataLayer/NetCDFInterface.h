#ifndef NETCDFINTERFACE_H
#define NETCDFINTERFACE_H

#include <stdio.h>
#include <netcdf.h>
#include "PhaseSpace.h"

#define ERR(e) {printf("NetCDF error: %s\n", nc_strerror(e));}


class NetCDFInterface
{
	public:
		// Constructor opens the file
		NetCDFInterface(const char* filename, FileMode mode);
		
		// create a new file and set up all the variables needed to write to it
		void openNewFileForWrite(const char* filename);
		void openFileForRead(const char* filename);
		
		// read a group of sequences
		void readSequencesAtCoordinate(Coord4 pos, Coord4 size);
		
		// write a group of sequences
		void writeSequence(Sequence& seq);
		
		// Destructor closes the file
		~NetCDFInterface();


	private:

		static const int NDIMS = 6;
		int dimids[NDIMS];
		int m_seqVarID;
		int m_ncid;
		FileMode m_mode;
		Count4D m_recordCounts;
	
};

#endif
