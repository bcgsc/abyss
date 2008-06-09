#ifndef IFILEREADER_H
#define IFILEREADER_H

#include <stdio.h>
#include <fstream>
#include "CommonDefs.h"
#include "Sequence.h"

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
