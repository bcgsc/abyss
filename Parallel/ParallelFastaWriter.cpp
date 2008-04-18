#include "ParallelFastaWriter.h"

ParallelFastaWriter::ParallelFastaWriter(const char* filename)
{	
	// Why doesnt mpi_open take in a const char*? const_cast anyway....
	MPI_File_open (MPI_COMM_WORLD, const_cast<char*>(filename), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &m_fileHandle);	
	MPI_File_set_view(m_fileHandle, 0, MPI_CHAR, MPI_CHAR, "native", MPI_INFO_NULL);
}

ParallelFastaWriter::~ParallelFastaWriter()
{	
	MPI_File_close(&m_fileHandle);
}

// Write out a single sequence
void ParallelFastaWriter::WriteSequence(const Sequence& seq, const int64_t id, const double multiplicity)
{
	int seqLength = seq.length();
	const int maxHeaderLength = 1000;
	int maxLength = seqLength + maxHeaderLength;
	printf("SEQ: %s\n", seq.c_str());
	char buffer[maxLength];
	int l = sprintf(buffer, ">%d %d %lf\n%s\n", (int)id, seqLength, multiplicity, seq.c_str());
	printf("BUFFER: %s\n", buffer);
	// Check for buffer overflow
	assert(l < maxLength);
	
	MPI_Status status;
	MPI_File_write_shared(m_fileHandle, buffer, l, MPI_CHAR, &status);
}
