#ifndef PAIRUTILS_H
#define PAIRUTILS_H 1

#include "Dictionary.h"

typedef std::string ContigID;
typedef unsigned LinearNumKey;

extern Dictionary g_contigIDs;

static inline LinearNumKey convertContigIDToLinearNumKey(
		const ContigID& id)
{
	return g_contigIDs.serial(id);
}

#endif
