/** Convert a graph from adj format to dot format.
 * Written by Shaun Jackman <sjackman@bcgsc.ca>.
 */
#include "ContigGraph.h"
#include <fstream>
#include <iostream>

using namespace std;

namespace opt {
	unsigned k; // used by Graph
}

int main(int argc, const char** argv)
{
	assert(argc <= 2);
	string path = argc > 1 ? argv[1] : "-";

	ifstream fin(path.c_str());
	if (path != "-")
		assert(fin.is_open());
	istream& in = path == "-" ? cin : fin;

	SimpleContigGraph g;
	in >> g;
	assert(in.eof());

	cout << "digraph \"" << path << "\" {\n"
		<< g
		<< "}\n";
	assert(cout.good());
	return 0;
}
