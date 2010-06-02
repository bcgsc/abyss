/** Convert a graph from adj format to dot format.
 * Written by Shaun Jackman <sjackman@bcgsc.ca>.
 */
#include "ContigGraph.h"
#include <iostream>

using namespace std;

namespace opt {
	unsigned k; // used by Graph
}

int main(int argc, const char** argv)
{
	assert(argc == 2);
	SimpleContigGraph g;
	loadGraphFromAdjFile(&g, argv[1]);
	cout << "digraph \"" << argv[1] << "\" {\n"
		<< g
		<< "}\n";
	assert(cout.good());
	return 0;
}
