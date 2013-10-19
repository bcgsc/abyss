#include "Graph/Path.h"
#include "Graph/BidirectionalBFS.h"
#include "Graph/BidirectionalBFSVisitor.h"
#include "Common/UnorderedMap.h"
#include "Common/Warnings.h"
#include <boost/graph/adjacency_list.hpp>
#include <gtest/gtest.h>

using namespace std;

typedef boost::adjacency_list<boost::vecS, boost::vecS,
	boost::bidirectionalS> Graph;

// note: vertex_descriptor for adjacency_list<> is int
typedef boost::graph_traits<Graph>::vertex_descriptor Vertex;
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

	void examine_vertex(Vertex u, const Graph& g, Direction dir)
	{
		SUPPRESS_UNUSED_WARNING(g);

#if 0
		// for debugging
		cerr << "visiting vertex: " << u << ", dir: " << dir << "\n";
#endif

		dirMap[u] = dir;
		rankMap[u] = rank++;
	}

	void common_edge(Edge e, const Graph& g)
	{
		SUPPRESS_UNUSED_WARNING(g);
		commonEdges.push_back(e);
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

	ASSERT_TRUE(commonEdges.size() == 1);
	ASSERT_TRUE(source(commonEdges[0], linearGraph) == 1u);
	ASSERT_TRUE(target(commonEdges[0], linearGraph) == 2u);
}

}
