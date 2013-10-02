/**
 * Connect pairs using a Bloom filter de Bruijn graph
 */

#include "DBGBloom.h"
#include "DBGBloomAlgorithms.h"

#include "Graph/DotIO.h"
#include "Graph/Options.h"

#include <cassert>

using namespace std;

namespace opt {
	unsigned k = 32;
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
	assert(argc > 1);
	DBGBloom g(opt::k);
	loadBloomFilter(g, argv[1]);

	write_dot(cout, g);

	processReads();

	return 0;
}
