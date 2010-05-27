#ifndef CONTIGID_H
#define CONTIGID_H 1

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
