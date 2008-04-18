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
		void WriteSequence(const Sequence& seq, const int64_t id, const double multiplicity);

	private:
		std::ofstream m_fileHandle;
};

#endif //FASTAWRITER_H
