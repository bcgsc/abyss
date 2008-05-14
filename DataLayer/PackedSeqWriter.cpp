#include "PackedSeqWriter.h"

PackedSeqWriter::PackedSeqWriter(const char* filename)
{	
	m_fileHandle.open(filename, std::ios::out | std::ios::binary);
	//printf("opening %s\n", filename);
	assert(m_fileHandle.is_open());		
}

PackedSeqWriter::~PackedSeqWriter()
{	
	m_fileHandle.close();
	assert(!m_fileHandle.is_open());
}

// Write out a single sequence
void PackedSeqWriter::WriteSequence(const PackedSeq& pSeq)
{
	// make sure the file is readable
	assert(m_fileHandle.is_open());
	
	const int size = sizeof(pSeq);
	char buffer[size];
	
	int nBytes = pSeq.serialize(buffer);

	assert(nBytes == size);
	(void)nBytes;
	// write as raw memory
	m_fileHandle.write(buffer, sizeof(pSeq));
}

void PackedSeqWriter::WriteSequence(const Sequence& pSeq)
{
	PackedSeq s(pSeq);
	WriteSequence(s);
}
