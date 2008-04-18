#include "FastaWriter.h"

FastaWriter::FastaWriter(const char* filename)
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
void FastaWriter::WriteSequence(const Sequence& seq, const int64_t id, const double multiplicity)
{
	
	// make sure the file is readable
	assert(m_fileHandle.is_open());

	m_fileHandle << ">" << id << " " << seq.length() << " " << multiplicity << "\n" << seq << "\n";

}
