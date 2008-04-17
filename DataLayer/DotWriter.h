#ifndef DOTWRITER_H
#define DOTWRITER_H 1

#include "ISequenceCollection.h"
#include <ostream>

class DotWriter {
	public:
		static void write(std::ostream& out,
				const ISequenceCollection& c);
};

#endif
