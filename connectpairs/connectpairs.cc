#include "DBGBloom.h"
#include "DBGBloomAlgorithms.h"
#include <cassert>

void buildDBGBloom() { assert(false); abort(); }

void findPaths() { assert(false); abort(); }

void mergePathts() { assert(false); abort(); }

void buildRead() { assert(false); abort(); }

void processRead()
{
	findPaths();
	mergePathts();
	buildRead();
}

void processReads()
{
	// for each read
	processRead();
}

int main()
{
	buildDBGBloom();
	processReads();

	return 0;
}
