#ifndef CONTIGGRAPH_H
#define CONTIGGRAPH_H 1

#include "DirectedGraph.h"

struct SimpleContigData
{
	unsigned length;
	SimpleContigData(unsigned length) : length(length) { }
	operator unsigned() const { return length; }
};

typedef DirectedGraph<SimpleContigData> SimpleContigGraph;

void loadGraphFromAdjFile(SimpleContigGraph* pGraph,
		const std::string& adjFile);

#endif
