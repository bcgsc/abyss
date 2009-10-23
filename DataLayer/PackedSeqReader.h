#ifndef PACKEDSEQREADER_H
#define PACKEDSEQREADER_H 1

#include "PackedSeq.h"
#include <fstream>
#include <vector>

class PackedSeqReader
{
	public:
		PackedSeqReader(const char* filename);
		~PackedSeqReader();
		bool ReadSequences(std::vector<PackedSeq>& outseqs);

	private:
		std::ifstream m_fileHandle;
		static const int m_numToRead = 1;
		int m_elementSize;
		int m_readSize;
		char* m_pBuffer;
};

#endif
