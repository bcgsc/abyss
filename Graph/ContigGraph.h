#ifndef CONTIGGRAPH_H
#define CONTIGGRAPH_H 1

#include "DirectedGraph.h"
#include <istream>

struct NoContigData { };

typedef DirectedGraph<NoContigData> ContigGraph;

std::istream& operator>>(std::istream& in, ContigGraph& o);

void readContigGraph(ContigGraph& graph, const std::string& path);

#endif
