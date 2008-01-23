#include "FastaWriter.h"

FastaWriter::FastaWriter(const char* filename) : m_count(0)
{	
	m_fileHandle.open(filename, std::ios::out);
	assert(m_fileHandle.is_open());
}

FastaWriter::~FastaWriter()
{	
	m_fileHandle.close();
	assert(!m_fileHandle.is_open());
}

// Write out a single sequence
void FastaWriter::WriteSequence(Sequence& seq)
{
	
	// make sure the file is readable
	assert(m_fileHandle.is_open());

	m_fileHandle << ">" << m_count << " " << seq.length() << "\n" << seq << "\n";
	
	m_count++;

}

// Write out a single sequence
void FastaWriter::WriteSequence(const PackedSeq& pSeq)
{
	
	// make sure the file is readable
	assert(m_fileHandle.is_open());

	m_fileHandle << ">" << m_count << " " << pSeq.getSequenceLength() << "\n" << pSeq.decode() << "\n";
	
	m_count++;

}

