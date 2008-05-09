#include "ParallelPackedSeqReader.h"

ParallelPackedSeqReader::ParallelPackedSeqReader(const char* filename, int /*id*/, int /*numNodes*/)
{
	
	//const int readUnit = sizeof(PackedSeq);
	printf("opening\n");
	
	// Why doesnt mpi_open take in a const char*? const_cast anyway....
	MPI_File_open (MPI_COMM_WORLD, const_cast<char*>(filename), MPI_MODE_RDONLY, MPI_INFO_NULL, &m_fileHandle);
	printf("opened\n");
	// Create the datatype we will read in
	// It is simply a byte stream with the size of the packedseq class
	// Todo: Make this network portable.
	const int numBytes = sizeof(PackedSeq);
	MPI_Type_contiguous(numBytes, MPI_BYTE, &m_readDatatype);
	
	//MPI_File_set_view(m_fileHandle, 0, m_readDatatype, m_readDatatype, "native", MPI_INFO_NULL);
}

ParallelPackedSeqReader::~ParallelPackedSeqReader()
{	
	MPI_File_close(&m_fileHandle);
}

int numRead = 0;
// Write out a single sequence
bool ParallelPackedSeqReader::ReadSequences(std::vector<Sequence>& seqs)
{
	MPI_Status status;
	const int numBytes = sizeof(PackedSeq);
	int count = 1000;
	char buffer[count*numBytes];
	MPI_File_read_shared(m_fileHandle, buffer, count, m_readDatatype, &status);
	
	int numActual;
	MPI_Get_count(&status, m_readDatatype, &numActual);
		
	for(int i = 0; i < numActual; i++)
	{
		PackedSeq s;
		s.unserialize((i*numBytes) + buffer);
		seqs.push_back(s.decode());
		
		numRead++;
		//printf("Read %s\n", s.decode().c_str());
	}
	
	if(numActual < count)
	{
		return false;
	}
	else
	{
		return true;
	}
}
