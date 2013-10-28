#include "Graph/Path.h"
#include "Graph/AllPathsSearch.h"
#include "Common/UnorderedMap.h"
#include "Common/Warnings.h"
#include <boost/graph/adjacency_list.hpp>
#include <gtest/gtest.h>
#include <set>

using namespace std;

typedef boost::adjacency_list<boost::vecS, boost::vecS,
	boost::bidirectionalS> Graph;

// note: vertex_descriptor for adjacency_list<> is int
typedef boost::graph_traits<Graph>::vertex_descriptor Vertex;
// note: edge_descriptor for adjacency_list<> is int
typedef boost::graph_traits<Graph>::edge_descriptor Edge;

namespace {

class AllPathsSearchTest : public ::testing::Test {

protected:

	Graph disconnectedGraph;
	Graph simpleAcyclicGraph;
	Graph cyclicGraph;

	AllPathsSearchTest() {

		add_edge(0, 1, disconnectedGraph);
		add_vertex(disconnectedGraph);

		add_edge(0, 1, simpleAcyclicGraph);
		add_edge(0, 2, simpleAcyclicGraph);
		add_edge(2, 3, simpleAcyclicGraph);

		add_edge(0, 1, cyclicGraph);
		add_edge(1, 2, cyclicGraph);
		add_edge(1, 3, cyclicGraph);
		add_edge(2, 3, cyclicGraph);
		add_edge(3, 4, cyclicGraph);
		add_edge(3, 5, cyclicGraph);
		add_edge(4, 5, cyclicGraph);
		add_edge(5, 6, cyclicGraph);

	}

};

TEST_F(AllPathsSearchTest, UnreachableGoal)
{
	vector< Path<Vertex> > paths;
	PathSearchResult result = allPathsSearch(disconnectedGraph, 0, 2, paths);

	EXPECT_EQ(NO_PATH, result);
	EXPECT_TRUE(paths.empty());
}

TEST_F(AllPathsSearchTest, StartNodeEqualsGoal)
{
	vector< Path<Vertex> > paths;
	PathSearchResult result = allPathsSearch(simpleAcyclicGraph, 0, 0, paths);

	EXPECT_EQ(FOUND_PATH, result);
	ASSERT_EQ(1u, paths.size());
	ASSERT_EQ("0", paths[0].str());
}

TEST_F(AllPathsSearchTest, SinglePath)
{
	vector< Path<Vertex> > paths;
	PathSearchResult result = allPathsSearch(simpleAcyclicGraph, 0, 3, paths);

	EXPECT_EQ(FOUND_PATH, result);
	ASSERT_EQ(1u, paths.size());
	ASSERT_EQ("0,2,3", paths[0].str());
}

TEST_F(AllPathsSearchTest, CyclicGraph)
{
	vector< Path<Vertex> > paths;
	PathSearchResult result = allPathsSearch(cyclicGraph, 0, 6, paths);

	set<string> expectedPaths;
	expectedPaths.insert("0,1,3,5,6");
	expectedPaths.insert("0,1,2,3,5,6");
	expectedPaths.insert("0,1,3,4,5,6");
	expectedPaths.insert("0,1,2,3,4,5,6");

	EXPECT_EQ(FOUND_PATH, result);
	ASSERT_EQ(4u, paths.size());

	// check that each path is one of the expected ones
	for (unsigned i = 0; i < 4; i++)
		ASSERT_TRUE(expectedPaths.find(paths[i].str()) != expectedPaths.end());

	// check that each path is unique
	for (unsigned i = 0; i < 3; i++)
		ASSERT_TRUE(paths[i].str() != paths[i+1].str());
}

}
