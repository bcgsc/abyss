#include "BloomDBG/bloom-dbg.h"
#include "BloomDBG/MaskedKmer.h"
#include "BloomDBG/RollingHash.h"

#include <gtest/gtest.h>
#include <iostream>

using namespace std;

TEST(BloomDBG, pathToSeq)
{
	const string inputSeq = "ACGTAC";
	const string spacedSeed = "10001";
	const unsigned k = 5;
	const unsigned numHashes = 2;

	MaskedKmer::setLength(k);
	MaskedKmer::setMask(spacedSeed);

	Path<BloomDBG::Vertex> path =
		BloomDBG::seqToPath(inputSeq, k, numHashes);
	ASSERT_EQ(2U, path.size());

	string outputSeq = BloomDBG::pathToSeq(path, k);
	ASSERT_EQ("ACNNAC", outputSeq);
}
