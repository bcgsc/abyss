#include "ContigGraph.h"
#include "ContigLength.h"
#include <cassert>
#include <cerrno>
#include <cstring> // for strerror
#include <fstream>
#include <limits> // for numeric_limits
#include <sstream>
#include <string>

using namespace std;

/** The length of each contig. */
static vector<unsigned> g_contigLengths;

/** Return the length of this contig in k-mer. */
unsigned ContigNode::length() const
{
	assert(!ambiguous());
	return g_contigLengths[id()];
}

static void assert_open(ifstream& f, const string& p)
{
	if (f.is_open())
		return;
	cerr << p << ": " << strerror(errno) << endl;
	exit(EXIT_FAILURE);
}

/** Read a contig graph. */
void readContigGraph(ContigGraph<>& graph, const std::string& path)
{
	ifstream fin(path.c_str());
	istream& in = path.empty() || path == "-" ? cin : fin;
	if (&in == &fin)
		assert_open(fin, path);
	assert(in.good());

	g_contigLengths = readContigLengths(in);
	in.clear();
	in.seekg(std::ios_base::beg);
	in >> graph;
	assert(in.eof());
}
