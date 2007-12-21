#ifndef PACKEDSEQWRITER_H
#define PACKEDSEQWRITER_H

#include <stdio.h>
#include <fstream>
#include "CommonDefs.h"
#include "Sequence.h"
#include "PackedSeq.h"

class PackedSeqWriter
{
	public:
	
		// Constructor opens file
		PackedSeqWriter(const char* filename, int sequenceLength);
		
		// Destructor closes it
		~PackedSeqWriter();
		
		// Write a single sequence
		void WriteSequence(PackedSeq& pSeq);
		
	private:
		std::ofstream m_fileHandle;
};

#endif //PACKEDSEQWRITER_H
