#ifndef PACKEDSEQWRITER_H
#define PACKEDSEQWRITER_H

#include "PackedSeq.h"
#include "Sequence.h"
#include <fstream>

class PackedSeqWriter
{
	public:
	
		// Constructor opens file
		PackedSeqWriter(const char* filename);
		
		// Destructor closes it
		~PackedSeqWriter();
		
		// Write a single sequence
		void WriteSequence(const PackedSeq& pSeq);
		void WriteSequence(const Sequence& pSeq);
		
	private:
		std::ofstream m_fileHandle;
};

#endif //PACKEDSEQWRITER_H
