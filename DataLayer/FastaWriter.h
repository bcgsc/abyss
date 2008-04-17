#ifndef FASTAWRITER_H
#define FASTAWRITER_H

#include <stdio.h>
#include <fstream>
#include "CommonDefs.h"
#include "Sequence.h"
#include "IFileWriter.h"
#include "PackedSeq.h"

class FastaWriter : public IFileWriter
{
	public:
	
		// Constructor opens file
		FastaWriter(const char* filename);
		
		// Destructor closes it
		~FastaWriter();
		
		// write a single sequence
		void WriteSequence(Sequence& seq, int64_t id, double multiplicity);
		void WriteSequence(const PackedSeq& pSeq, int64_t);
	private:
		std::ofstream m_fileHandle;
		int m_count;
};

#endif //FASTAWRITER_H
