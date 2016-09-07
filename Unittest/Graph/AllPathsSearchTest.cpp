#include "Graph/Path.h"
#include "Graph/AllPathsSearch.h"
#include "Common/UnorderedMap.h"
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
	Graph simpleCyclicGraph;
	Graph multiPathGraph;

	AllPathsSearchTest() {

		add_edge(0, 1, disconnectedGraph);
		add_vertex(disconnectedGraph);

		add_edge(0, 1, simpleAcyclicGraph);
		add_edge(0, 2, simpleAcyclicGraph);
		add_edge(2, 3, simpleAcyclicGraph);

		add_edge(0, 1, simpleCyclicGraph);
		add_edge(0, 4, simpleCyclicGraph);
		add_edge(1, 2, simpleCyclicGraph);
		add_edge(2, 1, simpleCyclicGraph);
		add_edge(1, 3, simpleCyclicGraph);

		add_edge(0, 1, multiPathGraph);
		add_edge(1, 2, multiPathGraph);
		add_edge(1, 3, multiPathGraph);
		add_edge(2, 3, multiPathGraph);
		add_edge(3, 4, multiPathGraph);
		add_edge(3, 5, multiPathGraph);
		add_edge(4, 5, multiPathGraph);
		add_edge(5, 6, multiPathGraph);

	}

};

TEST_F(AllPathsSearchTest, UnreachableGoal)
{
	AllPathsSearchResult<Vertex> result = allPathsSearch(disconnectedGraph, 0, 2);

	EXPECT_EQ(NO_PATH, result.resultCode);
	EXPECT_TRUE(result.paths.empty());
}

TEST_F(AllPathsSearchTest, StartNodeEqualsGoal)
{
	AllPathsSearchResult<Vertex> result = allPathsSearch(simpleAcyclicGraph, 0, 0);

	EXPECT_EQ(FOUND_PATH, result.resultCode);
	ASSERT_EQ(1u, result.paths.size());
	ASSERT_EQ("0", result.paths[0].str());
}

TEST_F(AllPathsSearchTest, SinglePath)
{
	AllPathsSearchResult<Vertex> result = allPathsSearch(simpleAcyclicGraph, 0, 3, 1, 2, 2, NO_LIMIT);

	EXPECT_EQ(FOUND_PATH, result.resultCode);
	ASSERT_EQ(1u, result.paths.size());
	ASSERT_EQ("0,2,3", result.paths[0].str());
}

TEST_F(AllPathsSearchTest, MultiPathGraph)
{
	AllPathsSearchResult<Vertex> result = allPathsSearch(multiPathGraph, 0, 6, 4, 4, 6, NO_LIMIT);

	set<string> expectedPaths;
	expectedPaths.insert("0,1,3,5,6");
	expectedPaths.insert("0,1,2,3,5,6");
	expectedPaths.insert("0,1,3,4,5,6");
	expectedPaths.insert("0,1,2,3,4,5,6");

	EXPECT_EQ(FOUND_PATH, result.resultCode);
	ASSERT_EQ(4u, result.paths.size());

	// check that each path is one of the expected ones
	for (unsigned i = 0; i < 4; i++)
		ASSERT_TRUE(expectedPaths.find(result.paths[i].str()) !=
			expectedPaths.end());

	// check that each path is unique
	for (unsigned i = 0; i < 3; i++)
		ASSERT_TRUE(result.paths[i].str() != result.paths[i+1].str());
}

TEST_F(AllPathsSearchTest, RespectsMaxPathsLimit)
{
	AllPathsSearchResult<Vertex> result = allPathsSearch(multiPathGraph, 0, 6, 3, NO_LIMIT, NO_LIMIT, NO_LIMIT);
	EXPECT_EQ(TOO_MANY_PATHS, result.resultCode);
}

TEST_F(AllPathsSearchTest, RespectsMaxDepthLimit)
{
	AllPathsSearchResult<Vertex> result = allPathsSearch(multiPathGraph, 0, 6, 4, 4, 5, NO_LIMIT);

	// We expect the fourth path ("0,1,2,3,4,5,6")
	// to be excluded by the max depth limit.  Note that
	// the depth of the start node is 0, and so a
	// path of length 7 reaches depth 6.

	set<string> expectedPaths;
	expectedPaths.insert("0,1,3,5,6");
	expectedPaths.insert("0,1,2,3,5,6");
	expectedPaths.insert("0,1,3,4,5,6");

	EXPECT_EQ(FOUND_PATH, result.resultCode);
	ASSERT_EQ(3u, result.paths.size());

	// check that each path is one of the expected ones
	for (unsigned i = 0; i < 3; i++)
		ASSERT_TRUE(expectedPaths.find(result.paths[i].str()) != expectedPaths.end());

	// check that each path is unique
	for (unsigned i = 0; i < 2; i++)
		ASSERT_TRUE(result.paths[i].str() != result.paths[i+1].str());
}

TEST_F(AllPathsSearchTest, RespectsMinDepthLimit)
{
	AllPathsSearchResult<Vertex> result = allPathsSearch(multiPathGraph, 0, 6, 4, 5, 6, NO_LIMIT);

	// We expect the shortest path ("0,1,3,4,6")
	// to be excluded by the min depth limit.  Note that
	// the depth of the start node is 0, and so a
	// path of length 5 reaches depth 4.

	set<string> expectedPaths;
	expectedPaths.insert("0,1,2,3,5,6");
	expectedPaths.insert("0,1,3,4,5,6");
	expectedPaths.insert("0,1,2,3,4,5,6");

	EXPECT_EQ(FOUND_PATH, result.resultCode);
	ASSERT_EQ(3u, result.paths.size());

	// check that each path is one of the expected ones
	for (unsigned i = 0; i < 3; i++)
		ASSERT_TRUE(expectedPaths.find(result.paths[i].str()) != expectedPaths.end());

	// check that each path is unique
	for (unsigned i = 0; i < 2; i++)
		ASSERT_TRUE(result.paths[i].str() != result.paths[i+1].str());
}

TEST_F(AllPathsSearchTest, PathContainsCycle)
{
	AllPathsSearchResult<Vertex> result = allPathsSearch(simpleCyclicGraph, 0, 3,
		NO_LIMIT, 0, NO_LIMIT, NO_LIMIT);
	EXPECT_EQ(PATH_CONTAINS_CYCLE, result.resultCode);
}

TEST_F(AllPathsSearchTest, IgnoreCycleNotOnPath)
{
	AllPathsSearchResult<Vertex> result = allPathsSearch(simpleCyclicGraph, 0, 4,
		NO_LIMIT, 0, NO_LIMIT, NO_LIMIT);

	EXPECT_EQ(FOUND_PATH, result.resultCode);
	ASSERT_EQ(1u, result.paths.size());
	ASSERT_EQ("0,4", result.paths.front().str());
}

}
