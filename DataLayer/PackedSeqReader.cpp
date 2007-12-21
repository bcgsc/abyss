#include "PackedSeqReader.h"

PackedSeqReader::PackedSeqReader(const char* filename)
{	
	m_fileHandle.open(filename, std::ios::binary);
	assert(m_fileHandle.is_open());
	
	// Read in the length of the sequences (the first record in the file)
	m_fileHandle >> m_seqLength; 
}

PackedSeqReader::~PackedSeqReader()
{	
	m_fileHandle.close();
	assert(!m_fileHandle.is_open());
}

bool PackedSeqReader::ReadAllSequences(PSequenceVector& outVector)
{
	// open file and check that we can read from it
	while(isGood())
	{
		PackedSeq* pSeq = ReadSequence();
		outVector.push_back(pSeq);
	}
}

// Read in a single sequence; this function allocates memory
PackedSeq* PackedSeqReader::ReadSequence()
{
	assert(m_fileHandle.is_open());
	
	// How many bytes need to be allocated
	int numBytes = PackedSeq::getNumCodingBytes(m_seqLength);
	
	// allocate storage for the data
	char* data = new char[numBytes];
	
	// read in the sequence
	m_fileHandle.read(data, numBytes);
		
	// allocate the new packed sequence
	PackedSeq* pSeq = new PackedSeq(data, m_seqLength);
	
	return pSeq;
}

bool PackedSeqReader::isGood()
{
	return !(m_fileHandle.eof() || m_fileHandle.peek() == EOF);
}
