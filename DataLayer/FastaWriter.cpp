#include "FastaWriter.h"
#include "Options.h"
#include <cassert>
#include <cstdio>
#include <cstdlib>

FastaWriter::FastaWriter(const char* path, bool append)
	: m_path(path), m_fileHandle(fopen(path, append ? "a" : "w"))
{
	if (m_fileHandle == NULL) {
		perror(m_path);
		exit(EXIT_FAILURE);
	}
}

FastaWriter::~FastaWriter()
{
	int n = fclose(m_fileHandle);
	if (n < 0) {
		perror(m_path);
		exit(EXIT_FAILURE);
	}
	m_fileHandle = NULL;
}

void FastaWriter::WriteSequence(const Sequence& seq,
		const int64_t id, const double multiplicity,
		const std::string& comment)
{
	assert(m_fileHandle != NULL);
	const char *sep = comment.empty() ? "" : " ";
	int n = opt::rank < 0
		? fprintf(m_fileHandle, ">%llu %zu %g%s%s\n%s\n",
				(long long unsigned)id,
				seq.length(), multiplicity,
				sep, comment.c_str(),
				seq.c_str())
		: fprintf(m_fileHandle, ">%u:%llu %zu %g%s%s\n%s\n",
				opt::rank,
				(long long unsigned)id,
				seq.length(), multiplicity,
				sep, comment.c_str(),
				seq.c_str());
	if (n < 0) {
		perror(m_path);
		exit(EXIT_FAILURE);
	}
}
