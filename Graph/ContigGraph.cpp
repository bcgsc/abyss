#include "ContigGraph.h"
#include "DirectedGraphImpl.h"
#include "PairUtils.h"
#include <cassert>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>

using namespace std;

// Explicit instantiation.
template class DirectedGraph<SimpleContigData>;

void parseAdjacencyLine(const string& adjLine, LinearNumKey currVert,
		SimpleContigGraph* pGraph);

/** Load an adjacency graph. */
void loadGraphFromAdjFile(SimpleContigGraph* pGraph,
		const string& lengthFile, const string& adjFile)
{
	// Load the lengths temporarily
	ContigLengthVec* pLengthVec = new ContigLengthVec();
	loadContigLengths(lengthFile, *pLengthVec);

	// First, load the vertices
	ifstream inStream(adjFile.c_str());
	assert(inStream.is_open());

	int numAdded = 0;
	LinearNumKey id;
	string adjRecord;
	while (inStream >> id
			&& getline(inStream, adjRecord)) {
		SimpleContigData data;
		data.length = pLengthVec->at(id);
		pGraph->addVertex(id, data);

		numAdded++;
		if (numAdded % 1000000 == 0)
			printf("added %d verts\n", numAdded);
	}
	assert(inStream.eof());

	// Delete the lengths to free up space
	delete pLengthVec;
	pLengthVec = NULL;

	// Now, load the edges
	inStream.clear();
	inStream.seekg(ios_base::beg);
	numAdded = 0;
	while (inStream >> id
			&& getline(inStream, adjRecord)) {
		parseAdjacencyLine(adjRecord, id, pGraph);

		numAdded++;
		if (numAdded % 1000000 == 0)
			printf("added edges for %d verts\n", numAdded);
	}
	assert(inStream.eof());

	size_t numVert = pGraph->getNumVertices();
	size_t numEdges = pGraph->countEdges(); // SLOW
	printf("Initial graph stats: num vert: %zu num edges: %zu\n",
			numVert, numEdges);
}

void parseAdjacencyLine(const string& adjLine, LinearNumKey currVert,
		SimpleContigGraph* pGraph)
{
	// convert to string stream
	stringstream ss(adjLine);
	for(size_t dirIdx = 0; dirIdx <= 1; ++dirIdx)
	{
		ss.ignore(numeric_limits<streamsize>::max(), '[');
		assert(ss.gcount() > 0);

		bool done = false;
		while(!done)
		{
			// extract the record
			string record;
			ss >> record;

			// check if the record is valid or we've hit the end
			if(record == "]")
			{
				done = true;
			} else {
				stringstream recSS(record);
				SimpleEdgeDesc sed;
				recSS >> sed;
				pGraph->addEdge(currVert,
						convertContigIDToLinearNumKey(sed.contig),
						(extDirection)dirIdx, sed.isRC);
			}
		}
	}
}
