#include "ContigGraph.h"
#include "ContigLength.h"
#include <cassert>
#include <fstream>
#include <limits> // for numeric_limits
#include <sstream>
#include <string>

using namespace std;

namespace opt {
	extern unsigned k;
};

/** The length of each contig. */
static vector<unsigned> g_contigLengths;

/** Return the length of this contig in k-mer. */
unsigned ContigNode::length() const
{
	assert(!ambiguous());
	return g_contigLengths[id()];
}

static void readEdges(istream& in, LinearNumKey id,
		ContigGraph& graph)
{
	for (int sense = false; sense <= true; ++sense) {
		string s;
		getline(in, s, !sense ? ';' : '\n');
		assert(in.good());
		istringstream ss(s);
		for (ContigNode edge; ss >> edge;)
			graph.add_edge(ContigNode(id, sense),
					sense ? ~edge : edge);
		assert(ss.eof());
	}
}

/** Read an adjacency graph. */
istream& operator>>(istream& in, ContigGraph& o)
{
	assert(g_contigLengths.empty());
	g_contigLengths = readContigLengths(in);

	// Create the vertices.
	o.clear();
	ContigGraph(2 * g_contigLengths.size()).swap(o);

	// Load the edges.
	in.clear();
	in.seekg(ios_base::beg);
	assert(in);
	for (string id; in >> id;) {
		in.ignore(numeric_limits<streamsize>::max(), ';');
		readEdges(in, stringToID(id), o);
	}
	assert(in.eof());
	return in;
}

/** Read an adjacency graph. */
void loadGraphFromAdjFile(ContigGraph* pGraph, const string& adjFile)
{
	ifstream in(adjFile.c_str());
	assert(in.is_open());
	in >> *pGraph;
	assert(in.eof());
}
