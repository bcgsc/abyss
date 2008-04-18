#ifndef PARALLELFASTAWRITER_H
#define PARALLELFASTAWRITER_H

#include <stdio.h>
#include <fstream>
#include <mpi.h>
#include "CommonDefs.h"
#include "Sequence.h"
#include "IFileWriter.h"
#include "PackedSeq.h"

class ParallelFastaWriter : public IFileWriter
{
	public:
	
		// Constructor opens file
		ParallelFastaWriter(const char* filename);
		
		// Destructor closes it
		~ParallelFastaWriter();
		
		// write a single sequence
		void WriteSequence(const Sequence& seq, const int64_t id, const double multiplicity);
		
	private:
		MPI_File m_fileHandle;
};

#endif //PARALLELFASTAWRITER_H
