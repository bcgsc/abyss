#include "PackedSeqWriter.h"

PackedSeqWriter::PackedSeqWriter(const char* filename, int sequenceLength)
{	
	m_fileHandle.open(filename, std::ios::out | std::ios::binary);
	//printf("opening %s\n", filename);
	assert(m_fileHandle.is_open());	
	// write out the length of the sequences
	m_fileHandle.write((char*)&sequenceLength, sizeof(int));
	
}

PackedSeqWriter::~PackedSeqWriter()
{	
	m_fileHandle.close();
	assert(!m_fileHandle.is_open());
}

// Write out a single sequence
void PackedSeqWriter::WriteSequence(const PackedSeq& pSeq, int64_t id)
{
	
	// make sure the file is readable
	assert(m_fileHandle.is_open());
	
	m_fileHandle.write(pSeq.getDataPtr(), pSeq.getNumCodingBytes(pSeq.getSequenceLength()));
}
