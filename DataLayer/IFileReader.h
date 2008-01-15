#ifndef IFILEREADER_H
#define IFILEREADER_H

#include <stdio.h>
#include <fstream>
#include "CommonDefs.h"
#include "Sequence.h"
#include "PackedSeq.h"

class IFileReader
{
	public:
		// Read in a single sequence
		virtual ~IFileReader() {};
		virtual PackedSeq ReadSequence() = 0;
		virtual bool isGood() = 0;
};

#endif //IFILEREADER_H
