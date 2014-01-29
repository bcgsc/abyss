#ifndef CONNECTPAIRS_H
#define CONNECTPAIRS_H
#include "DataLayer/FastaInterleave.h"

static unsigned maskNew(const FastqRecord& read1, const FastqRecord& read2,
		FastaRecord& merged);

#endif
