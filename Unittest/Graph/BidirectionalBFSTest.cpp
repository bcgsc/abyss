#include "Graph/Path.h"
#include "Graph/BidirectionalBFS.h"
#include "Graph/BidirectionalBFSVisitor.h"
#include "Common/UnorderedMap.h"
#include "Common/Warnings.h"
#include <boost/graph/adjacency_list.hpp>
#include <gtest/gtest.h>

using namespace std;

static const bool DEBUG = false;

typedef boost::adjacency_list<boost::vecS, boost::vecS,
	boost::bidirectionalS> Graph;

// note: vertex_descriptor for adjacency_list<> is int
typedef boost::graph_traits<Graph>::vertex_descriptor Vertex;
// note: edge_descriptor for adjacency_list<> is int
typedef boost::graph_traits<Graph>::edge_descriptor Edge;

class TestVisitor : public BidirectionalBFSVisitor<Graph>
{
public:

	typedef unordered_map<Vertex, Direction> DirMap;
	typedef unordered_map<Vertex, int> RankMap;
	typedef vector<Edge> EdgeList;

	DirMap dirMap;
	RankMap rankMap;
	EdgeList commonEdges;
	int rank;

	TestVisitor() : rank(0) { }

	BFSVisitorResult discover_vertex(const Vertex& u, const Graph& g, Direction dir,
		unsigned numFrontierNodes)
	{
		if (DEBUG)
			cerr << dir << ": discover_vertex " << u << "\n";

		SUPPRESS_UNUSED_WARNING(g);
		SUPPRESS_UNUSED_WARNING(numFrontierNodes);
		return SUCCESS;
	}

	void examine_vertex(const Vertex& u, const Graph& g, Direction dir)
	{
		if (DEBUG)
			cerr << dir << ": examine_vertex " << u << "\n";

		SUPPRESS_UNUSED_WARNING(g);
		dirMap[u] = dir;
		rankMap[u] = rank++;
	}

	void examine_edge(const Edge& e, const Graph& g, Direction dir)
	{
		if (DEBUG)
			cerr << dir << ": examine_edge " << e << "\n";

		SUPPRESS_UNUSED_WARNING(g);
	}

	BFSVisitorResult tree_edge(const Edge& e, const Graph& g, Direction dir)
	{
		if (DEBUG)
			cerr << dir << ": tree_edge " << e << "\n";

		SUPPRESS_UNUSED_WARNING(g);
		return SUCCESS;
	}

	BFSVisitorResult common_edge(const Edge& e, const Graph& g, Direction dir)
	{
		if (DEBUG)
			cerr << dir << ": common_edge " << e << "\n";

		SUPPRESS_UNUSED_WARNING(g);
		commonEdges.push_back(e);
		return SUCCESS;
	}

	BFSVisitorResult non_tree_edge(const Edge& e, const Graph& g, Direction dir)
	{
		if (DEBUG)
			cerr << dir << ": non_tree_edge " << e << "\n";

		SUPPRESS_UNUSED_WARNING(g);
		commonEdges.push_back(e);
		return SUCCESS;
	}

	void gray_target(const Edge& e, const Graph& g, Direction dir)
	{
		if (DEBUG)
			cerr << dir << ": gray_target " << e << "\n";

		SUPPRESS_UNUSED_WARNING(g);
	}

	void black_target(const Edge& e, const Graph& g, Direction dir)
	{
		if (DEBUG)
			cerr << dir << ": black_target " << e << "\n";

		SUPPRESS_UNUSED_WARNING(g);
	}

};

namespace {

class BidirectionalBFSTest : public ::testing::Test {

protected:

	Graph linearGraph;
	Graph branchingGraph;

	BidirectionalBFSTest() {

		add_edge(0, 1, linearGraph);
		add_edge(1, 2, linearGraph);
		add_edge(2, 3, linearGraph);

		add_edge(0, 1, branchingGraph);
		add_edge(0, 2, branchingGraph);
		add_edge(1, 3, branchingGraph);
		add_edge(3, 4, branchingGraph);
		add_edge(4, 6, branchingGraph);
		add_edge(5, 6, branchingGraph);

	}

};

TEST_F(BidirectionalBFSTest, AlternatesDirection)
{
	TestVisitor visitor;
	bidirectionalBFS(linearGraph, 0, 3, visitor);

	TestVisitor::DirMap& dir = visitor.dirMap;

	ASSERT_EQ(dir[0], FORWARD);
	ASSERT_EQ(dir[3], REVERSE);
	ASSERT_EQ(dir[1], FORWARD);
	ASSERT_EQ(dir[2], REVERSE);
}

TEST_F(BidirectionalBFSTest, FollowsBreadthFirstOrder)
{
	TestVisitor visitor;
	bidirectionalBFS(branchingGraph, 0, 6, visitor);

	TestVisitor::RankMap& rank = visitor.rankMap;

	ASSERT_TRUE(rank[1] > rank[0]);
	ASSERT_TRUE(rank[2] > rank[0]);
	ASSERT_TRUE(rank[3] > rank[1]);
	ASSERT_TRUE(rank[3] > rank[2]);
	ASSERT_TRUE(rank[3] > rank[4]);
	ASSERT_TRUE(rank[3] > rank[5]);
	ASSERT_TRUE(rank[4] > rank[6]);
	ASSERT_TRUE(rank[5] > rank[6]);
}

TEST_F(BidirectionalBFSTest, IdentifiesCommonEdge)
{
	TestVisitor visitor;
	bidirectionalBFS(linearGraph, 0, 3, visitor);

	TestVisitor::EdgeList& commonEdges = visitor.commonEdges;

	// Note: Each common edge is visited twice.  It is
	// visited once by the forward traversal and once
	// by the reverse traversal.

	ASSERT_TRUE(commonEdges.size() == 2);
	ASSERT_TRUE(commonEdges[0] == commonEdges[1]);
	ASSERT_TRUE(source(commonEdges[0], linearGraph) == 1u);
	ASSERT_TRUE(target(commonEdges[0], linearGraph) == 2u);
}

}
