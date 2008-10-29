#include "FastaWriter.h"
#include "Options.h"
#include <cstdio>

FastaWriter::FastaWriter(const char* filename, bool append)
	: m_fileHandle(fopen(filename, append ? "a" : "w"))
{
	assert(m_fileHandle != NULL);
}

FastaWriter::~FastaWriter()
{
	fclose(m_fileHandle);
	m_fileHandle = NULL;
}

// Write out a single sequence
void FastaWriter::WriteSequence(const Sequence& seq,
		const int64_t id, const double multiplicity)
{
	assert(m_fileHandle != NULL);
	if (opt::rank < 0) {
		fprintf(m_fileHandle, ">%llu %zu %g\n%s\n",
				(long long unsigned)id,
				seq.length(), multiplicity, seq.c_str());
	} else {
		fprintf(m_fileHandle, ">%u:%llu %zu %g\n%s\n", opt::rank,
				(long long unsigned)id,
				seq.length(), multiplicity, seq.c_str());
	}
}
