#ifndef PARALLELFASTAREADER_H
#define PARALLELFASTAREADER_H

#include <stdio.h>
#include <fstream>
#include <mpi.h>
#include "CommonDefs.h"
#include "Sequence.h"
#include "IFileWriter.h"
#include "PackedSeq.h"

class ParallelPackedSeqReader
{
	public:
	
		// Constructor opens file
		ParallelPackedSeqReader(const char* filename, int id, int numNodes);
		
		// Destructor closes it
		~ParallelPackedSeqReader();
		
		// read a single sequence
		bool ReadSequences(std::vector<Sequence>& seqs);
		
	private:
		MPI_File m_fileHandle;
		MPI_Datatype m_readDatatype;
};

#endif //PARALLELFASTAWRITER_H
