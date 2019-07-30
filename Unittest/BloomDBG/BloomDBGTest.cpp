#include "Common/Sequence.h"
#include "BloomDBG/bloom-dbg.h"
#include "BloomDBG/MaskedKmer.h"
#include "BloomDBG/RollingHash.h"
#include "BloomDBG/RollingBloomDBG.h"
#include "vendor/btl_bloomfilter/BloomFilter.hpp"

#include <gtest/gtest.h>
#include <iostream>

using namespace std;
typedef RollingBloomDBG<BloomFilter> Graph;
typedef graph_traits<Graph> GraphTraits;

/* each vertex is represented by
 * std::pair<MaskedKmer, vector<size_t>>, where 'string' is the
 * k-mer and 'vector<size_t>' is the associated set of
 * hash values */
typedef graph_traits<Graph>::vertex_descriptor V;

/** Convert a path in the de Bruijn graph to a sequence */
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
