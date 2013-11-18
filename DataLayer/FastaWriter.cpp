#include "FastaWriter.h"
#include "Common/Options.h"
#include <cassert>
#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring> // for strerror
#include <iostream>
#include <unistd.h> // for fsync

using namespace std;

static inline void die(const string& s)
{
	cerr << "error: writing to `" << s << "': "
		<< strerror(errno) << endl;
	exit(EXIT_FAILURE);
}

FastaWriter::FastaWriter(const char* path, bool append)
	: m_path(path), m_fileHandle(fopen(path, append ? "a" : "w"))
{
	if (m_fileHandle == NULL)
		die(m_path);
}

FastaWriter::~FastaWriter()
{
	int n = fsync(fileno(m_fileHandle));
	if (n < 0)
		die(m_path);
	n = fclose(m_fileHandle);
	if (n < 0)
		die(m_path);
	m_fileHandle = NULL;
}

void FastaWriter::WriteSequence(const Sequence& seq, unsigned id,
		unsigned multiplicity, const string& comment)
{
	assert(m_fileHandle != NULL);
	const char *sep = comment.empty() ? "" : " ";
	int n = opt::rank < 0
		? fprintf(m_fileHandle, ">%llu %zu %u%s%s\n%s\n",
				(long long unsigned)id,
				seq.length(), multiplicity,
				sep, comment.c_str(),
				seq.c_str())
		: fprintf(m_fileHandle, ">%u:%llu %zu %u%s%s\n%s\n",
				opt::rank,
				(long long unsigned)id,
				seq.length(), multiplicity,
				sep, comment.c_str(),
				seq.c_str());
	if (n < 0)
		die(m_path);
}

void FastaWriter::WriteSequence(const Sequence& seq, unsigned long long id, const std::string& comment)
{
	assert(m_fileHandle != NULL);
	int n = fprintf(m_fileHandle, ">%llu %s\n%s\n", id, comment.c_str(), seq.c_str());
	if (n < 0)
		die(m_path);
}

void FastaWriter::WriteSequence(const Sequence& seq, const std::string& id, const std::string& comment)
{
	assert(m_fileHandle != NULL);
	int n = fprintf(m_fileHandle, ">%s %s\n%s\n", id.c_str(), comment.c_str(), seq.c_str());
	if (n < 0)
		die(m_path);
}
