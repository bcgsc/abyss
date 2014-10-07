#ifndef PAIREDDBG_DOTWRITER_H
#define PAIREDDBG_DOTWRITER_H 1

#include "SequenceCollection.h"
#include <ostream>

class DotWriter {
  public:
	static void write(std::ostream& out,
			const SequenceCollectionHash& c);
};

#endif
