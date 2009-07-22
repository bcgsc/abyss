#ifndef CONTIGGRAPH_H
#define CONTIGGRAPH_H 1

#include "DirectedGraph.h"

struct SimpleContigData
{
	int length;
};

typedef DirectedGraph<SimpleContigData> SimpleContigGraph;

void loadGraphFromAdjFile(SimpleContigGraph* pGraph,
		std::string& lengthFile, std::string adjFile);

void parseAdjacencyLine(std::string& adjLine, LinearNumKey currVert,
		SimpleContigGraph* pGraph);

#endif
