#include "FastaWriter.h"
#include "Options.h"
#include <cstdio>
#include <unistd.h>

FastaWriter::FastaWriter(const char* filename, bool append)
	: m_fileHandle(fopen(filename, append ? "a" : "w"))
{
	assert(m_fileHandle != NULL);
	if (!append)
		freopen(NULL, "a", m_fileHandle);
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
		// Non-parallel
		fprintf(m_fileHandle, ">%llu %u %g\n%s\n",
				id, seq.length(), multiplicity, seq.c_str());
	} else {
		// Parallel
		int fd = fileno(m_fileHandle);
		int f_lock = lockf(fd, F_LOCK, 0);
		assert(f_lock == 0);
		(void)f_lock;
		fprintf(m_fileHandle, ">%u:%llu %u %g\n%s\n", opt::rank,
				id, seq.length(), multiplicity, seq.c_str());
		fflush(m_fileHandle);
		lseek(fd, 0, SEEK_SET);
		int f_ulock = lockf(fd, F_ULOCK, 0);
		assert(f_ulock == 0);
		(void)f_ulock;
	}
}
