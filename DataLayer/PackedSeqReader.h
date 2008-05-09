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
		virtual bool ReadSequences(PSequenceVector& outseqs);
		
		
	private:
		std::ifstream m_fileHandle;
		static const int m_numToRead = 1;
		int m_elementSize;
		int m_readSize;
		char* m_pBuffer;
};

#endif //PACKEDSEQREADER_H
