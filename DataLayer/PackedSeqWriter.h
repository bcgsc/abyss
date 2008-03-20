#ifndef PACKEDSEQWRITER_H
#define PACKEDSEQWRITER_H

#include <stdio.h>
#include <fstream>
#include "CommonDefs.h"
#include "Sequence.h"
#include "IFileWriter.h"
#include "PackedSeq.h"

class PackedSeqWriter : public IFileWriter
{
	public:
	
		// Constructor opens file
		PackedSeqWriter(const char* filename, int sequenceLength);
		
		// Destructor closes it
		~PackedSeqWriter();
		
		// Write a single sequence
		void WriteSequence(const PackedSeq& pSeq, int64_t id = 0);
		
	private:
		std::ofstream m_fileHandle;
};

#endif //PACKEDSEQWRITER_H
