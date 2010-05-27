#include "ContigGraph.h"
#include "ContigID.h"
#include "DirectedGraphImpl.h"
#include <cassert>
#include <fstream>
#include <iostream>
#include <limits> // for numeric_limits
#include <sstream>
#include <string>

using namespace std;

namespace opt {
	/** Abort the search after visiting maxCost vertices. */
	unsigned maxCost = 100000;
};

// Explicit instantiation.
template class DirectedGraph<SimpleContigData>;

static void readEdges(istream& in, LinearNumKey id,
		SimpleContigGraph& graph)
{
	for (extDirection dir = SENSE; dir <= ANTISENSE; ++dir) {
		string s;
		getline(in, s, dir == SENSE ? ';' : '\n');
		assert(in.good());
		istringstream ss(s);
		for (ContigNode edge; ss >> edge;)
			graph.addEdge(id, dir, edge);
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
		assert(length >= opt::k);
		pGraph->addVertex(stringToID(id), length);
		if (++count % 1000000 == 0)
			cout << "Read " << count << " vertices" << endl;
	}
	assert(in.eof());

	// Load the edges.
	in.clear();
	in.seekg(ios_base::beg);
	count = 0;
	while (in >> id) {
		in.ignore(numeric_limits<streamsize>::max(), ';');
		readEdges(in, stringToID(id), *pGraph);
		if (++count % 1000000 == 0)
			cout << "Read edges for " << count << " vertices" << endl;
	}
	assert(in.eof());
}
