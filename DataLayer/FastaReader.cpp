#include "FastaReader.h"

FastaReader::FastaReader(const char* filename)
{	
	m_fileHandle.open(filename);
	assert(m_fileHandle.is_open());
}

FastaReader::~FastaReader()
{	
	m_fileHandle.close();
	assert(!m_fileHandle.is_open());
}

bool FastaReader::ReadAllSequences(PSequenceVector& outVector)
{
	// open file and check that we can read from it
	while(isGood())
	{
		PackedSeq* pSeq = ReadSequence();	
		outVector.push_back(pSeq);
	}
}

// Read in a single sequence to the out parameter, return whether there are more sequences to read
PackedSeq* FastaReader::ReadSequence()
{
	char headerBuffer[MAX_FASTA_LINE];
	char seqBuffer[MAX_FASTA_LINE];	
	char id[SEQUENCE_ID_LENGTH];
	
	// make sure the file is readable
	assert(m_fileHandle.is_open());

	
	// read in the header
	m_fileHandle.getline(headerBuffer, MAX_FASTA_LINE);

	// read in the sequence
	m_fileHandle.getline(seqBuffer, MAX_FASTA_LINE);
		
	// check if the reads were successful
	if(m_fileHandle.fail())
	{
		printf("error reading from fasta file. Is the line length too long?\n");
		return false;
	}
	
	// parse the header
	if(sscanf(headerBuffer, ">%s %*s", id) != 1)
	{
		printf("invalid header format, read failed\n");
		assert(false);
	}
	
	Sequence seq(seqBuffer);
	
	// Create the packed seq
	PackedSeq* pSeq = new PackedSeq(seq);
	return pSeq;

}

bool FastaReader::isGood()
{
	return !(m_fileHandle.eof() || m_fileHandle.peek() == EOF);
}
