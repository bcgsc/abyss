#ifndef FASTAWRITER_H
#define FASTAWRITER_H

#include "IFileWriter.h"
#include "Sequence.h"
#include <cstdio>

class FastaWriter : public IFileWriter
{
	public:
	
		// Constructor opens file
		FastaWriter(const char* filename, bool append = false);
		
		// Destructor closes it
		~FastaWriter();
		
		// write a single sequence
		void WriteSequence(const Sequence& seq, const int64_t id, const double multiplicity);

	private:
		FILE* m_fileHandle;
};

#endif //FASTAWRITER_H
