#ifndef IFILEWRITER_H
#define IFILEWRITER_H

#include <stdio.h>
#include <fstream>
#include "CommonDefs.h"
#include "Sequence.h"
#include "PackedSeq.h"

class IFileWriter
{
	public:
		// Read in a single sequence
		virtual ~IFileWriter() {};
		virtual void WriteSequence(const PackedSeq& pSeq, int64_t id) = 0;
};

#endif //IFILEWRITER_H
