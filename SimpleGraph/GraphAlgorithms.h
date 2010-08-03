#ifndef GRAPHALGORITHMS_H
#define GRAPHALGORITHMS_H 1

#include "ContigGraph.h"
#include "ContigNode.h"
#include "ContigPath.h"
#include <istream>
#include <utility>
#include <vector>

namespace opt {
	extern unsigned k;
	extern unsigned maxCost;

	/** Abort the search after visiting maxPaths solutions. */
	static const unsigned maxPaths = 200;
}

/** Contig properties. */
struct Contig {
	unsigned length;

	friend std::istream& operator >>(std::istream& in, Contig& o)
	{
		if (in >> o.length) {
			assert(o.length >= opt::k);
			o.length -= opt::k - 1;
		}
		return in;
	}
};

typedef ContigGraph<Contig> Graph;
typedef std::pair<ContigNode, unsigned> Constraint;
typedef std::vector<Constraint> Constraints;
typedef std::vector<ContigPath> ContigPaths;

bool depthFirstSearch(const Graph& g,
		ContigNode origin, Constraints& constraints,
		ContigPaths& superPaths, unsigned& compCost);

#endif
