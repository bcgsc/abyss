#ifndef CONTIGGRAPH_H
#define CONTIGGRAPH_H 1

#include "DirectedGraph.h"
#include <istream>

struct SimpleContigData { };

typedef DirectedGraph<SimpleContigData> SimpleContigGraph;

std::istream& operator>>(std::istream& in, SimpleContigGraph& o);

void loadGraphFromAdjFile(SimpleContigGraph* pGraph,
		const std::string& adjFile);

#endif
