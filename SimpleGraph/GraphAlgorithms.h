#ifndef GRAPHALGORITHMS_H
#define GRAPHALGORITHMS_H 1

#include "ContigGraph.h"
#include "ContigNode.h"
#include "ContigPath.h"
#include <utility>
#include <vector>

namespace opt {
	extern unsigned k;
	extern unsigned maxCost;

	/** Abort the search after visiting maxPaths solutions. */
	static const unsigned maxPaths = 200;
}

typedef std::pair<ContigNode, unsigned> Constraint;
typedef std::vector<Constraint> Constraints;
typedef std::vector<ContigPath> ContigPaths;

bool depthFirstSearch(const ContigGraph<>& g,
		const ContigNode& sourceKey, Constraints& constraints,
		ContigPaths& superPaths, unsigned& compCost);

#endif
