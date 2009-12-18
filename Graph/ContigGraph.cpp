#include "ContigGraph.h"
#include "DirectedGraphImpl.h"
#include "PairUtils.h"
#include <cassert>
#include <fstream>
#include <iostream>
#include <limits> // for numeric_limits
#include <sstream>
#include <string>

using namespace std;

// Explicit instantiation.
template class DirectedGraph<SimpleContigData>;

static void readEdges(istream& in, LinearNumKey id,
		SimpleContigGraph& graph)
{
	for (extDirection dir = SENSE; dir <= ANTISENSE; ++dir) {
		string s;
		getline(in, s, dir == SENSE ? ';' : '\n');
		cout << s << endl;
		assert(in.good());
		istringstream ss(s);
		for (SimpleEdgeDesc edge; ss >> edge;)
			graph.addEdge(id,
					convertContigIDToLinearNumKey(edge.contig),
					dir, edge.isRC);
		assert(ss.eof());
	}
}

/** Load an adjacency graph. */
void loadGraphFromAdjFile(SimpleContigGraph* pGraph,
		const string& adjFile)
{
	// Load the vertices.
	ifstream in(adjFile.c_str());
	assert(in.is_open());

	unsigned count = 0;
	string id;
	unsigned length;
	while (in >> id >> length) {
		in.ignore(numeric_limits<streamsize>::max(), '\n');
		pGraph->addVertex(convertContigIDToLinearNumKey(id), length);
		if (++count % 1000000 == 0)
			cout << "Read " << count << " vertices" << endl;
	}
	assert(in.eof());

	// Load the edges.
	in.clear();
	in.seekg(ios_base::beg);
	count = 0;
	while (in >> id >> length) {
		readEdges(in, convertContigIDToLinearNumKey(id), *pGraph);
		if (++count % 1000000 == 0)
			cout << "Read edges for " << count << " vertices" << endl;
	}
	assert(in.eof());
}
