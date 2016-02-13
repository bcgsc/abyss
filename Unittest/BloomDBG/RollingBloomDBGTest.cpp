#include "BloomDBG/RollingBloomDBG.h"
#include "BloomDBG/RollingHashIterator.h"
#include "lib/bloomfilter-2dfba08d120d7659e8c75cf5c501b3b9040e98cb/BloomFilter.hpp"
#include "Common/UnorderedSet.h"

#include <gtest/gtest.h>
#include <string>

using namespace std;
using namespace boost;

typedef RollingBloomDBG<BloomFilter> Graph;
typedef graph_traits<Graph> GraphTraits;

/* each vertex is represented by
 * std::pair<Kmer, vector<size_t>>, where 'string' is the
 * k-mer and 'vector<size_t>' is the associated set of
 * hash values */
typedef graph_traits<Graph>::vertex_descriptor V;

/** Test fixture for RollingBloomDBG tests. */
class RollingBloomDBGTest : public ::testing::Test
{
protected:

	const unsigned m_k;
	const unsigned m_bloomSize;
	const unsigned m_bloomHashes;
	BloomFilter m_bloom;
	Graph m_graph;

	RollingBloomDBGTest() : m_k(5), m_bloomSize(100000), m_bloomHashes(2),
		m_bloom(m_bloomSize, m_bloomHashes, m_k), m_graph(m_bloom)
	{
		Kmer::setLength(m_k);

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

		m_bloom.insert("CGACT");
		m_bloom.insert("TGACT");
		m_bloom.insert("GACTC");
		m_bloom.insert("ACTCT");
		m_bloom.insert("ACTCG");
	}

};

TEST_F(RollingBloomDBGTest, out_edge_iterator)
{
	/* TEST: check that "GACTC" has the expected outgoing edges */

	const V GACTC(Kmer("GACTC"), RollingHash("GACTC", m_bloomHashes, m_k));
	const V ACTCT(Kmer("ACTCT"), RollingHash("ACTCT", m_bloomHashes, m_k));
	const V ACTCG(Kmer("ACTCG"), RollingHash("ACTCG", m_bloomHashes, m_k));

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

	const V GACTC(Kmer("GACTC"), RollingHash("GACTC", m_bloomHashes, m_k));
	const V ACTCT(Kmer("ACTCT"), RollingHash("ACTCT", m_bloomHashes, m_k));
	const V ACTCG(Kmer("ACTCG"), RollingHash("ACTCG", m_bloomHashes, m_k));

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

	const V GACTC(Kmer("GACTC"), RollingHash("GACTC", m_bloomHashes, m_k));
	const V CGACT(Kmer("CGACT"), RollingHash("CGACT", m_bloomHashes, m_k));
	const V TGACT(Kmer("TGACT"), RollingHash("TGACT", m_bloomHashes, m_k));

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

/** Test fixture for RollingBloomDBG with spaced seed k-mers. */
class RollingBloomDBGSpacedSeedTest : public ::testing::Test
{
protected:

	const unsigned m_k;
	const unsigned m_bloomSize;
	const unsigned m_bloomHashes;
	BloomFilter m_bloom;
	Graph m_graph;
	const std::string m_spacedSeed;

	RollingBloomDBGSpacedSeedTest() : m_k(5), m_bloomSize(100000), m_bloomHashes(2),
		m_bloom(m_bloomSize, m_bloomHashes, m_k), m_graph(m_bloom),
		m_spacedSeed("11011")
	{
		Kmer::setLength(m_k);

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

		m_bloom.insert(RollingHash("CGACT", m_bloomHashes,
			m_k, m_spacedSeed).getHash());
		m_bloom.insert(RollingHash("TGACT", m_bloomHashes,
			m_k, m_spacedSeed).getHash());
		m_bloom.insert(RollingHash("GACTC", m_bloomHashes,
			m_k, m_spacedSeed).getHash());
		m_bloom.insert(RollingHash("ACTCT", m_bloomHashes,
			m_k, m_spacedSeed).getHash());
		m_bloom.insert(RollingHash("ACTCG", m_bloomHashes,
			m_k, m_spacedSeed).getHash());
	}

};

TEST_F(RollingBloomDBGSpacedSeedTest, out_edge_iterator)
{
	/* TEST: check that "GACTC" has the expected outgoing edges */

	const V GACTC(Kmer("GACTC"),
		RollingHash("GACTC", m_bloomHashes, m_k, m_spacedSeed));
	const V ACTCT(Kmer("ACTCT"),
		RollingHash("ACTCT", m_bloomHashes, m_k, m_spacedSeed));
	const V ACTCG(Kmer("ACTCG"),
		RollingHash("ACTCG", m_bloomHashes, m_k, m_spacedSeed));

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
