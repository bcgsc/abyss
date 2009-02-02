#ifndef IFILEWRITER_H
#define IFILEWRITER_H 1

#include "Sequence.h"

class IFileWriter
{
	public:
		// Read in a single sequence
		virtual ~IFileWriter() {};
		virtual void WriteSequence(const Sequence& seq, const int64_t id, const double multiplicity) = 0;
};

#endif //IFILEWRITER_H
