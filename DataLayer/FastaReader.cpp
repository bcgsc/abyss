#include "FastaReader.h"
#include "Options.h"

FastaReader::FastaReader(const char* filename)
{	
	m_fileHandle.open(filename);
	assert(m_fileHandle.is_open());
}

FastaReader::~FastaReader()
{	
	m_fileHandle.close();
	assert(!m_fileHandle.is_open());
}


Sequence FastaReader::ReadSequence()
{
	char headerBuffer[MAX_FASTA_LINE];
	char seqBuffer[MAX_FASTA_LINE];	
	char id[MAX_FASTA_LINE];
	
	// make sure the file is readable
	assert(m_fileHandle.is_open());

	// read in the header
	m_fileHandle.getline(headerBuffer, MAX_FASTA_LINE);

	// read in the sequence
	m_fileHandle.getline(seqBuffer, MAX_FASTA_LINE);
			
	// parse the header
	if(sscanf(headerBuffer, ">%s %*s", id) != 1)
	{
		printf("invalid header format, read failed\n");
		assert(false);
	}
	
	return Sequence(seqBuffer);	
}

// Read in a group of sequences and return whether there are sequences remaining
bool FastaReader::ReadSequences(PSequenceVector& outseqs)
{
	Sequence seq = ReadSequence();
	size_t pos = seq.find_first_not_of("ACGT");
	if (pos == std::string::npos) {
		outseqs.push_back(seq);
	} else {
		if (opt::verbose > 0)
			fprintf(stderr,
					"warning: discarded sequence containing `%c'\n",
					seq[pos]);
		m_nonacgt++;
	}
	return isGood();
}

bool FastaReader::isGood()
{
	return !(m_fileHandle.eof() || m_fileHandle.peek() == EOF);
}
