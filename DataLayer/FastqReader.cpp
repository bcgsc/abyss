#include "FastqReader.h"
#include "Options.h"
#include <algorithm>
#include <cctype>

FastqReader::FastqReader(const char* filename)
	: m_nonacgt(0)
{
	m_fileHandle.open(filename);
	assert(m_fileHandle.is_open());
}

FastqReader::~FastqReader()
{
	m_fileHandle.close();
}


Sequence FastqReader::ReadSequence()
{
	char buf[MAX_FASTA_LINE];
	char seq[MAX_FASTA_LINE];

	// Read the header.
	m_fileHandle.getline(buf, sizeof buf);
	assert(m_fileHandle.gcount() < MAX_FASTA_LINE-1);
	assert(buf[0] == '@');

	// Read the sequence.
	m_fileHandle.getline(seq, sizeof seq);
	assert(m_fileHandle.gcount() < MAX_FASTA_LINE-1);

	// Read the quality values.
	m_fileHandle.getline(buf, sizeof buf);
	assert(m_fileHandle.gcount() < MAX_FASTA_LINE-1);
	assert(buf[0] == '+');
	m_fileHandle.getline(buf, sizeof buf);
	assert(m_fileHandle.gcount() < MAX_FASTA_LINE-1);

	Sequence s(seq);
	transform(s.begin(), s.end(), s.begin(), ::toupper);
	return s;
}

/**
 * Read in a group of sequences and return whether any sequences
 * remain.
 */
bool FastqReader::ReadSequences(SequenceVector& outseqs)
{
	Sequence seq = ReadSequence();
	size_t pos = seq.find_first_not_of("ACGT");
	if (pos == std::string::npos) {
		outseqs.push_back(seq);
	} else {
		if (opt::verbose > 1)
			fprintf(stderr,
					"warning: discarded sequence containing `%c'\n",
					seq[pos]);
		m_nonacgt++;
	}
	return isGood();
}

bool FastqReader::isGood()
{
	return !(m_fileHandle.eof() || m_fileHandle.peek() == EOF);
}
