#ifndef PACKEDSEQREADER_H
#define PACKEDSEQREADER_H

#include "IFileReader.h"
#include "PackedSeq.h"
#include <fstream>

class PackedSeqReader /*: public IFileReader*/
{
	public:
	
		// Constructor opens file
		PackedSeqReader(const char* filename);
		
		// Destructor closes it
		~PackedSeqReader();
		
		// Read in a single sequence
		/*virtual*/ bool ReadSequences(PSequenceVector& outseqs);
		
		// Returns the number of sequences containing non-ACGT
		// characters, which is impossible for a packed sequence.
		/*virtual unsigned getNonACGT() { return 0; }*/

	private:
		std::ifstream m_fileHandle;
		static const int m_numToRead = 1;
		int m_elementSize;
		int m_readSize;
		char* m_pBuffer;
};

#endif //PACKEDSEQREADER_H
