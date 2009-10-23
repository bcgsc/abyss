#include "PackedSeqReader.h"
#include <cassert>
#include <cstdio>

using namespace std;

PackedSeqReader::PackedSeqReader(const char* filename)
{	
	m_fileHandle.open(filename, ios::binary);
	m_elementSize = PackedSeq::serialSize();
	m_readSize = m_numToRead * m_elementSize;
	m_pBuffer = new char[m_readSize];
		
	assert(m_fileHandle.is_open());
	if (m_fileHandle.peek() == EOF)
		fprintf(stderr, "warning: `%s' is empty\n", filename);
}

PackedSeqReader::~PackedSeqReader()
{	
	// free the buffer
	delete [] m_pBuffer;
	m_pBuffer = 0;
	
	// close the file
	m_fileHandle.close();
	assert(!m_fileHandle.is_open());
}

bool PackedSeqReader::ReadSequences(vector<PackedSeq>& outseqs)
{
	assert(m_fileHandle.is_open());
	if (m_fileHandle.eof())
		return false;

	// read in the sequence
	m_fileHandle.read(m_pBuffer, m_readSize);
	int numBytesRead = m_fileHandle.gcount();

	// Calculate the number of sequences read in
	int numSequencesRead = numBytesRead / m_elementSize;

	// ensure it is a whole number
	assert(numBytesRead % m_elementSize == 0);	
	
	// build the sequences
	for(int i = 0; i < numSequencesRead; i++)
	{
		PackedSeq seq;
		seq.unserialize(i*m_elementSize + m_pBuffer);
		outseqs.push_back(seq);
	}

	return true;
}
