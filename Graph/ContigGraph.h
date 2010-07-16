#ifndef CONTIGGRAPH_H
#define CONTIGGRAPH_H 1

#include "DirectedGraph.h"
#include <istream>

struct NoContigData { };

class ContigGraph;

std::istream& operator>>(std::istream& in, ContigGraph& o);

void readContigGraph(ContigGraph& graph, const std::string& path);

class ContigGraph : public DirectedGraph<NoContigData> {
	typedef DirectedGraph<NoContigData> DG;

  public:
	/** Construct an empty contig graph. */
	ContigGraph() { }

	/** Construct a contig graph with n vertices. */
	ContigGraph(vertices_size_type n) : DG(n) { }

  private:
	ContigGraph(const ContigGraph&);
};

#endif
