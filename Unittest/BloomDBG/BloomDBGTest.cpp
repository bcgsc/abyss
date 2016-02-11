#include "BloomDBG/bloom-dbg.h"
#include "Common/Kmer.h"
#include "BloomDBG/RollingHash.h"

#include <gtest/gtest.h>
#include <iostream>

using namespace std;

TEST(BloomDBG, pathToSeq)
{
	typedef std::pair<Kmer, RollingHash> V;

	const string inputSeq = "ACGTAC";
	const string spacedSeed = "10001";
	const unsigned k = 5;
	const unsigned numHashes = 2;

	Kmer::setLength(k);

	Path<V> path = BloomDBG::seqToPath(inputSeq, k, numHashes, spacedSeed);
	ASSERT_EQ(2U, path.size());

	string outputSeq = BloomDBG::pathToSeq(path, k, spacedSeed);
	ASSERT_EQ("ACNNAC", outputSeq);
}
