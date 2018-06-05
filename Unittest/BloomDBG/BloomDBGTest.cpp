#include "Common/Sequence.h"
#include "BloomDBG/bloom-dbg.h"
#include "BloomDBG/MaskedKmer.h"
#include "BloomDBG/RollingHash.h"
#include "BloomDBG/RollingBloomDBG.h"
#include "lib/bloomfilter/BloomFilter.hpp"

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

/** Split a sequence at branching k-mers the de Bruijn graph */
TEST(BloomDBG, splitSeq)
{
	const size_t bloomSize = 100000;
	const unsigned k = 5;
	const unsigned numHashes = 2;
	const unsigned minBranchLen = 1;
	size_t hashes[MAX_HASHES];

	/* it is important to reset these, since they persist between tests */
	MaskedKmer::setLength(k);
	MaskedKmer::mask().clear();

	/*
	 * Test graph (k=5):
	 *
	 *   GACTC-ACTCG-CTCGG
	 *
	 * Input sequence (horizontal path above):
	 *
	 *   GACTCGG
	 */

	BloomFilter bloom1(bloomSize, numHashes, k);

	RollingHash("GACTC", numHashes, k).getHashes(hashes);
	bloom1.insert(hashes);
	RollingHash("ACTCG", numHashes, k).getHashes(hashes);
	bloom1.insert(hashes);
	RollingHash("CTCGG", numHashes, k).getHashes(hashes);
	bloom1.insert(hashes);

	Sequence seq1 = "GACTCGG";

	Graph graph1(bloom1);
	vector<Sequence> segments1 = BloomDBG::splitSeq(seq1, k,
		numHashes, graph1, minBranchLen);

	V GACTC(V("GACTC", RollingHash("GACTC", numHashes, k)));

	ASSERT_EQ(1U, out_degree(GACTC, graph1));
	ASSERT_EQ(1U, segments1.size());
	ASSERT_EQ("GACTCGG", segments1.front());

	/*
	 * Test graph (k=5):
	 *
	 *         ACTCT
	 *        /
	 *   GACTC-ACTCG-CTCGG
	 *              /
	 *         TCTCG
	 *
	 * Input sequence (horizontal path above):
	 *
	 *   GACTCGG
	 */

	BloomFilter bloom2(bloomSize, numHashes, k);

	RollingHash("GACTC", numHashes, k).getHashes(hashes);
	bloom2.insert(hashes);
	RollingHash("ACTCT", numHashes, k).getHashes(hashes);
	bloom2.insert(hashes);
	RollingHash("ACTCG", numHashes, k).getHashes(hashes);
	bloom2.insert(hashes);
	RollingHash("CTCGG", numHashes, k).getHashes(hashes);
	bloom2.insert(hashes);
	RollingHash("TCTCG", numHashes, k).getHashes(hashes);
	bloom2.insert(hashes);

	Sequence seq2 = "GACTCGG";

	Graph graph2(bloom2);
	vector<Sequence> segments2 = BloomDBG::splitSeq(seq2, k,
		numHashes, graph2, minBranchLen);

	ASSERT_EQ(3U, segments2.size());
	ASSERT_EQ("GACTC", segments2.at(0));
	ASSERT_EQ("GACTCGG", segments2.at(1));
	ASSERT_EQ("CTCGG", segments2.at(2));

	/*
	 * Test graph (k=5):
	 *
	 *   TACTC       CTCGA
	 *        \     /
	 *   GACTC-ACTCG-CTCGG
	 *
	 * Input sequence (horizontal path above):
	 *
	 *   ACTCG
	 */

	BloomFilter bloom3(bloomSize, numHashes, k);

	RollingHash("TACTC", numHashes, k).getHashes(hashes);
	bloom2.insert(hashes);
	RollingHash("GACTC", numHashes, k).getHashes(hashes);
	bloom2.insert(hashes);
	RollingHash("ACTCG", numHashes, k).getHashes(hashes);
	bloom2.insert(hashes);
	RollingHash("CTCGA", numHashes, k).getHashes(hashes);
	bloom2.insert(hashes);
	RollingHash("CTCGG", numHashes, k).getHashes(hashes);
	bloom2.insert(hashes);

	Sequence seq3 = "ACTCG";

	Graph graph3(bloom3);
	vector<Sequence> segments3 = BloomDBG::splitSeq(seq3, k,
		numHashes, graph3, minBranchLen);

	ASSERT_EQ(1U, segments3.size());
	ASSERT_EQ("ACTCG", segments3.front());
}
