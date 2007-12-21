#ifndef PACKEDSEQREADER_H
#define PACKEDSEQREADER_H

#include <stdio.h>
#include <fstream>
#include "IFileReader.h"
#include "CommonDefs.h"
#include "Sequence.h"
#include "PackedSeq.h"

class PackedSeqReader : public IFileReader
{
	public:
	
		// Constructor opens file
		PackedSeqReader(const char* filename);
		
		// Destructor closes it
		~PackedSeqReader();
		
		// Read in a single sequence
		PackedSeq* ReadSequence();
		
		// Read all sequences in the file
		bool ReadAllSequences(PSequenceVector& outVector);
		
		// Returns true if there are sequences left to read
		bool isGood();
		
	private:
		std::ifstream m_fileHandle;
		int m_seqLength;
};

#endif //PACKEDSEQREADER_H
