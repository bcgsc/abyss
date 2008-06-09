#ifndef FASTQREADER_H
#define FASTQREADER_H 1

#include "CommonDefs.h"
#include "IFileReader.h"
#include "Sequence.h"
#include <fstream>

class FastqReader : public IFileReader
{
	public:
		FastqReader(const char* filename);
		~FastqReader();
		Sequence ReadSequence();
		virtual bool ReadSequences(SequenceVector& outseqs);
		bool isGood();
		virtual unsigned getNonACGT() { return m_nonacgt; }

	private:
		std::ifstream m_fileHandle;
		unsigned m_nonacgt;
};

#endif
