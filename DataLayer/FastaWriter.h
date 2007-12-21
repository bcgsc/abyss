#ifndef FASTAWRITER_H
#define FASTAWRITER_H

#include <stdio.h>
#include <fstream>
#include "CommonDefs.h"
#include "Sequence.h"
#include "PackedSeq.h"

class FastaWriter
{
	public:
	
		// Constructor opens file
		FastaWriter(const char* filename);
		
		// Destructor closes it
		~FastaWriter();
		
		// Read in a single sequence to the out parameter, return whether there are more sequences to read
		void WriteSequence(Sequence& seq);
		
	private:
		std::ofstream m_fileHandle;
		int m_count;
};

#endif //FASTAWRITER_H
