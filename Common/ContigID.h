#ifndef CONTIGID_H
#define CONTIGID_H 1

#include "Dictionary.h"

typedef std::string ContigID;
typedef unsigned LinearNumKey;

extern Dictionary g_contigIDs;

/** Convert a string to a numeric contig ID. */
static inline LinearNumKey convertContigIDToLinearNumKey(
		const ContigID& id)
{
	return g_contigIDs.serial(id);
}

/** Convert a numeric contig ID to a string. */
static inline const std::string& idToString(unsigned id)
{
	return g_contigIDs.key(id);
}

#endif
