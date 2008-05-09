#ifndef FASTAREADER_H
#define FASTAREADER_H

#include <stdio.h>
#include <fstream>
#include "IFileReader.h"
#include "CommonDefs.h"
#include "Sequence.h"
#include "PackedSeq.h"

class FastaReader : public IFileReader
{
	public:
	
		// Constructor opens file
		FastaReader(const char* filename);
		
		// Destructor closes it
		~FastaReader();
		
		// Read a single sequence from the file
		Sequence ReadSequence();
		
		// Read sequences into the vector as packed seqs
		// Returns true unless eof has been reached
		virtual bool ReadSequences(PSequenceVector& outseqs);
		
		// Returns true unless eof has been reached
		bool isGood();
				
	private:

		std::ifstream m_fileHandle;
};

#endif //FASTAREADER_H
