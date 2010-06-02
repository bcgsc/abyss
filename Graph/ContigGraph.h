#ifndef CONTIGGRAPH_H
#define CONTIGGRAPH_H 1

#include "DirectedGraph.h"

struct SimpleContigData { };

typedef DirectedGraph<SimpleContigData> SimpleContigGraph;

void loadGraphFromAdjFile(SimpleContigGraph* pGraph,
		const std::string& adjFile);

#endif
