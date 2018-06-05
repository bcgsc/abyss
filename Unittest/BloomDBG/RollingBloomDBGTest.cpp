#include "BloomDBG/RollingBloomDBG.h"
#include "lib/bloomfilter/BloomFilter.hpp"
#include "Common/UnorderedSet.h"

#include <gtest/gtest.h>
#include <string>

using namespace std;
using namespace boost;

typedef RollingBloomDBG<BloomFilter> Graph;
typedef graph_traits<Graph> GraphTraits;
typedef graph_traits<Graph>::vertex_descriptor V;

/** Test fixture for RollingBloomDBG tests. */
class RollingBloomDBGTest : public ::testing::Test
{
protected:

	const unsigned m_k;
	const unsigned m_bloomSize;
	const unsigned m_numHashes;
	BloomFilter m_bloom;
	Graph m_graph;

	RollingBloomDBGTest() : m_k(5), m_bloomSize(100000), m_numHashes(2),
		m_bloom(m_bloomSize, m_numHashes, m_k), m_graph(m_bloom)
	{
		MaskedKmer::setLength(m_k);

		/*
		 * Test de Bruijn graph:
		 *
		 *  CGACT       ACTCT
		 *       \     /
		 *        GACTC
		 *       /     \
		 *  TGACT       ACTCG
		 *
		 * Note: No unexpected edges
		 * are created by the reverse
		 * complements of these k-mers.
		 */

		size_t hashes[MAX_HASHES];
		RollingHash("CGACT", m_numHashes, m_k).getHashes(hashes);
		m_bloom.insert(hashes);
		RollingHash("TGACT", m_numHashes, m_k).getHashes(hashes);
		m_bloom.insert(hashes);
		RollingHash("GACTC", m_numHashes, m_k).getHashes(hashes);
		m_bloom.insert(hashes);
		RollingHash("ACTCT", m_numHashes, m_k).getHashes(hashes);
		m_bloom.insert(hashes);
		RollingHash("ACTCG", m_numHashes, m_k).getHashes(hashes);
		m_bloom.insert(hashes);
	}

};

TEST_F(RollingBloomDBGTest, out_edge_iterator)
{
	/* TEST: check that "GACTC" has the expected outgoing edges */

	const V GACTC("GACTC", RollingHash("GACTC", m_numHashes, m_k));
	const V ACTCT("ACTCT", RollingHash("ACTCT", m_numHashes, m_k));
	const V ACTCG("ACTCG", RollingHash("ACTCG", m_numHashes, m_k));

	unordered_set<V> expectedNeighbours;
	expectedNeighbours.insert(ACTCT);
	expectedNeighbours.insert(ACTCG);

	ASSERT_EQ(2u, out_degree(GACTC, m_graph));
	GraphTraits::out_edge_iterator ei, ei_end;
	boost::tie(ei, ei_end) = out_edges(GACTC, m_graph);
	ASSERT_NE(ei_end, ei);
	unordered_set<V>::iterator neighbour =
		expectedNeighbours.find(target(*ei, m_graph));
	EXPECT_NE(expectedNeighbours.end(), neighbour);
	expectedNeighbours.erase(neighbour);
	ei++;
	ASSERT_NE(ei_end, ei);
	neighbour = expectedNeighbours.find(target(*ei, m_graph));
	ASSERT_NE(expectedNeighbours.end(), neighbour);
	ei++;
	ASSERT_EQ(ei_end, ei);
}

TEST_F(RollingBloomDBGTest, adjacency_iterator)
{
	/* TEST: check that "GACTC" has the expected outgoing edges */

	const V GACTC("GACTC", RollingHash("GACTC", m_numHashes, m_k));
	const V ACTCT("ACTCT", RollingHash("ACTCT", m_numHashes, m_k));
	const V ACTCG("ACTCG", RollingHash("ACTCG", m_numHashes, m_k));

	unordered_set<V> expectedNeighbours;
	expectedNeighbours.insert(ACTCT);
	expectedNeighbours.insert(ACTCG);

	ASSERT_EQ(2u, out_degree(GACTC, m_graph));
	GraphTraits::adjacency_iterator ai, ai_end;
	boost::tie(ai, ai_end) = adjacent_vertices(GACTC, m_graph);
	ASSERT_NE(ai_end, ai);
	unordered_set<V>::iterator neighbour =
		expectedNeighbours.find(*ai);
	EXPECT_NE(expectedNeighbours.end(), neighbour);
	expectedNeighbours.erase(neighbour);
	ai++;
	ASSERT_NE(ai_end, ai);
	neighbour = expectedNeighbours.find(*ai);
	ASSERT_NE(expectedNeighbours.end(), neighbour);
	ai++;
	ASSERT_EQ(ai_end, ai);
}

TEST_F(RollingBloomDBGTest, in_edges)
{
	/* TEST: check that "GACTC" has the expected ingoing edges */

	const V GACTC("GACTC", RollingHash("GACTC", m_numHashes, m_k));
	const V CGACT("CGACT", RollingHash("CGACT", m_numHashes, m_k));
	const V TGACT("TGACT", RollingHash("TGACT", m_numHashes, m_k));

	unordered_set<V> expectedNeighbours;
	expectedNeighbours.insert(CGACT);
	expectedNeighbours.insert(TGACT);

	ASSERT_EQ(2u, in_degree(GACTC, m_graph));
	GraphTraits::in_edge_iterator ei, ei_end;
	boost::tie(ei, ei_end) = in_edges(GACTC, m_graph);
	ASSERT_NE(ei_end, ei);
	unordered_set<V>::iterator neighbour =
		expectedNeighbours.find(source(*ei, m_graph));
	EXPECT_NE(expectedNeighbours.end(), neighbour);
	expectedNeighbours.erase(neighbour);
	ei++;
	ASSERT_NE(ei_end, ei);
	neighbour = expectedNeighbours.find(source(*ei, m_graph));
	ASSERT_NE(expectedNeighbours.end(), neighbour);
	ei++;
	ASSERT_EQ(ei_end, ei);
}

TEST_F(RollingBloomDBGTest, pathTraversal)
{
	/*
	 * Walk a simple path:
	 *
	 * CGACT-GACTC-ACTCG
	 */

	BloomFilter bloom(m_bloomSize, m_numHashes, m_k);
	Graph graph(bloom);

	const V CGACT("CGACT", RollingHash("CGACT", m_numHashes, m_k));
	const V GACTC("GACTC", RollingHash("GACTC", m_numHashes, m_k));
	const V ACTCG("ACTCG", RollingHash("ACTCG", m_numHashes, m_k));

	size_t hashes[MAX_HASHES];
	CGACT.rollingHash().getHashes(hashes);
	bloom.insert(hashes);
	GACTC.rollingHash().getHashes(hashes);
	bloom.insert(hashes);
	ACTCG.rollingHash().getHashes(hashes);
	bloom.insert(hashes);

	/* step one */

	V v = CGACT;
    ASSERT_EQ(1u, out_degree(v, graph));
	GraphTraits::out_edge_iterator ei, ei_end;
	boost::tie(ei, ei_end) = out_edges(v, graph);
	ASSERT_NE(ei_end, ei);
	ASSERT_EQ(CGACT, source(*ei, graph));
	ASSERT_EQ(GACTC, target(*ei, graph));
	v = target(*ei, graph);
	++ei;
	ASSERT_EQ(ei_end, ei);

	/* step two */

    ASSERT_EQ(1u, out_degree(v, graph));
	boost::tie(ei, ei_end) = out_edges(v, graph);
	ASSERT_NE(ei_end, ei);
	ASSERT_EQ(GACTC, source(*ei, graph));
	ASSERT_EQ(ACTCG, target(*ei, graph));
	v = target(*ei, graph);
	++ei;
	ASSERT_EQ(ei_end, ei);
}

/** Test fixture for RollingBloomDBG with spaced seed k-mers. */
class RollingBloomDBGSpacedSeedTest : public ::testing::Test
{
protected:

	const unsigned m_k;
	const unsigned m_bloomSize;
	const unsigned m_numHashes;
	BloomFilter m_bloom;
	Graph m_graph;
	const std::string m_spacedSeed;

	RollingBloomDBGSpacedSeedTest() : m_k(5), m_bloomSize(100000), m_numHashes(1),
		m_bloom(m_bloomSize, m_numHashes, m_k), m_graph(m_bloom),
		m_spacedSeed("11011")
	{
		MaskedKmer::setLength(m_k);
		MaskedKmer::setMask(m_spacedSeed);

		/*
		 * Test de Bruijn graph:
		 *
		 *  CGACT       ACTCT
		 *       \     /
		 *        GACTC
		 *       /     \
		 *  TGACT       ACTCG
		 *
		 * Masked version:
		 *
		 *  CG_CT       AC_CT
		 *       \     /
		 *        GA_TC
		 *       /     \
		 *  TG_CT       AC_CG
		 *
		 * Note: With respect to the spaced seed "11011",
		 * GACTC is equivalent to its own reverse complement
		 * GAGTC.  However, this does not result in
		 * any additional edges in the graph.
		 */

		size_t hashes[MAX_HASHES];
		RollingHash("CGACT", m_numHashes, m_k).getHashes(hashes);
		m_bloom.insert(hashes);
		RollingHash("TGACT", m_numHashes, m_k).getHashes(hashes);
		m_bloom.insert(hashes);
		RollingHash("GACTC", m_numHashes, m_k).getHashes(hashes);
		m_bloom.insert(hashes);
		RollingHash("ACTCT", m_numHashes, m_k).getHashes(hashes);
		m_bloom.insert(hashes);
		RollingHash("ACTCG", m_numHashes, m_k).getHashes(hashes);
		m_bloom.insert(hashes);
	}

};

TEST_F(RollingBloomDBGSpacedSeedTest, out_edge_iterator)
{
	/* TEST: check that "GACTC" has the expected outgoing edges */

	const V GACTC("GACTC", RollingHash("GACTC", m_numHashes, m_k));
	const V ACTCT("ACTCT", RollingHash("ACTCT", m_numHashes, m_k));
	const V ACTCG("ACTCG", RollingHash("ACTCG", m_numHashes, m_k));

	unordered_set<V> expectedNeighbours;
	expectedNeighbours.insert(ACTCT);
	expectedNeighbours.insert(ACTCG);

	ASSERT_EQ(2u, out_degree(GACTC, m_graph));
	GraphTraits::out_edge_iterator ei, ei_end;
	boost::tie(ei, ei_end) = out_edges(GACTC, m_graph);
	ASSERT_NE(ei_end, ei);
	unordered_set<V>::iterator neighbour =
		expectedNeighbours.find(target(*ei, m_graph));
	EXPECT_NE(expectedNeighbours.end(), neighbour);
	expectedNeighbours.erase(neighbour);
	ei++;
	ASSERT_NE(ei_end, ei);
	neighbour = expectedNeighbours.find(target(*ei, m_graph));
	ASSERT_NE(expectedNeighbours.end(), neighbour);
	ei++;
	ASSERT_EQ(ei_end, ei);
}
