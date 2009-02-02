#ifndef IFILEREADER_H
#define IFILEREADER_H

#include "Sequence.h"

const unsigned MAX_FASTA_LINE = 262144;

class IFileReader
{
	public:
		// Read in a single sequence
		virtual ~IFileReader() {};
		
		// Returns true if there are more sequences to read
		virtual bool ReadSequences(SequenceVector& outseqs) = 0;

		// Returns the number of sequences containing non-ACGT
		// characters.
		virtual unsigned getNonACGT() = 0;
};

#endif //IFILEREADER_H
