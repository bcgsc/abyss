/**
 * Connect pairs using a Bloom filter de Bruijn graph
 */

#include "DBGBloom.h"
#include "DBGBloomAlgorithms.h"

#include "Common/Options.h"
#include "Graph/DotIO.h"
#include "Graph/Options.h"
#include "Graph/GraphUtil.h"

#include <cassert>

using namespace std;

namespace opt {
	/** The size of a k-mer. */
	unsigned k = 16;
}

/** Load the bloom filter. */
static void loadBloomFilter(DBGBloom& g, const string& path)
{
	g.open(path);
}

static void findPaths() { assert(false); abort(); }

static void mergePathts() { assert(false); abort(); }

static void buildRead() { assert(false); abort(); }

static void processRead()
{
	findPaths();
	mergePathts();
	buildRead();
}

static void processReads()
{
	// for each read
	processRead();
}

/**
 * Connect pairs using a Bloom filter de Bruijn graph
 */
int main(int argc, const char* argv[])
{
	opt::verbose = 1;

	assert(opt::k > 0);
	Kmer::setLength(opt::k);

	assert(argc > 1);
	DBGBloom g(opt::k);
	loadBloomFilter(g, argv[1]);
	if (opt::verbose > 0)
		printGraphStats(cerr, g);

	write_dot(cout, g);

	processReads();

	return 0;
}
