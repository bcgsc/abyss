#include "FastaWriter.h"

FastaWriter::FastaWriter(const char* filename, bool append)
{	
	std::ios_base::openmode mode = std::ios::out;
	if(append)
	{
		// if the append flag is set, open for append
		mode |= std::ios::app;
	}
	else
	{
		// otherwise truncate any current data in the file
		mode |= std::ios::trunc;
	}
	
	m_fileHandle.open(filename, mode);
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
