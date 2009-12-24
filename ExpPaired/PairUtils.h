#ifndef PAIRUTILS_H
#define PAIRUTILS_H 1

#include "ContigNode.h"
#include "Dictionary.h"
#include <vector>

typedef std::string ContigID;
typedef unsigned LinearNumKey;
typedef ContigNode SimpleEdgeDesc;

extern Dictionary g_contigIDs;

static inline LinearNumKey convertContigIDToLinearNumKey(
		const ContigID& id)
{
	return g_contigIDs.serial(id);
}

typedef std::vector<int> ContigLengthVec;

void loadContigLengths(const std::string& path,
		ContigLengthVec& lengths);

#endif
