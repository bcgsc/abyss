#include "FastaReader.h"
#include "Options.h"
#include <algorithm>
#include <cassert>
#include <cctype>

FastaReader::FastaReader(const char* filename)
	: m_nonacgt(0)
{	
	m_fileHandle.open(filename);
	assert(m_fileHandle.is_open());
	if (m_fileHandle.peek() == EOF)
		fprintf(stderr, "warning: `%s' is empty\n", filename);
}

FastaReader::~FastaReader()
{	
	m_fileHandle.close();
	assert(!m_fileHandle.is_open());
}


Sequence FastaReader::ReadSequence()
{
	char buf[MAX_FASTA_LINE];

	// Discard comments.
	while (m_fileHandle.peek() == '#') {
		m_fileHandle.getline(buf, (ssize_t)sizeof buf);
		assert(m_fileHandle.gcount() < (ssize_t)sizeof buf - 1);
		if (m_fileHandle.peek() == EOF) {
			fputs("error: file ends in comments\n", stderr);
			assert(false);
		}
	}

	// Read the header.
	m_fileHandle.getline(buf, (ssize_t)sizeof buf);
	assert(m_fileHandle.gcount() < (ssize_t)sizeof buf - 1);
	assert(buf[0] == '>');

	// Read the sequence.
	m_fileHandle.getline(buf, (ssize_t)sizeof buf);
	assert(m_fileHandle.gcount() < (ssize_t)sizeof buf - 1);

	Sequence s(buf);
	transform(s.begin(), s.end(), s.begin(), ::toupper);
	return s;
}

// Read in a group of sequences and return whether there are sequences remaining
bool FastaReader::ReadSequences(SequenceVector& outseqs)
{
	if (!isGood())
		return false;
	Sequence seq = ReadSequence();
	size_t pos = seq.find_first_not_of("ACGT");
	if (pos == std::string::npos) {
		outseqs.push_back(seq);
	} else {
		if (opt::verbose > 3)
			fprintf(stderr,
					"warning: discarded sequence containing `%c'\n",
					seq[pos]);
		m_nonacgt++;
	}
	return true;
}

bool FastaReader::isGood()
{
	return !(m_fileHandle.eof() || m_fileHandle.peek() == EOF);
}
