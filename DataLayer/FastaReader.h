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
		
		// Read in a single sequence
		PackedSeq ReadSequence();
		
		// Read all sequences in the file
		bool ReadAllSequences(PSequenceVector& outVector);
		
		// Returns true if there are sequences left to read
		bool isGood();
		
	private:
		std::ifstream m_fileHandle;
};

#endif //FASTAREADER_H
