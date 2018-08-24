#include "Graph/Path.h"
#include "Graph/BreadthFirstSearch.h"
#include "Graph/ConstrainedBFSVisitor.h"
#include "Graph/DefaultColorMap.h"
#include <boost/graph/adjacency_list.hpp>
#include <boost/property_map/property_map.hpp>
#include <gtest/gtest.h>

using namespace std;
using namespace boost;

typedef adjacency_list<vecS, vecS, bidirectionalS> Graph;
typedef graph_traits<Graph>::vertex_descriptor V;
typedef std::vector< Path<V> > PathList;

namespace {

class ConstrainedBFSVisitorTest : public ::testing::Test {

protected:

	Graph simpleAcyclicGraph;
	Graph simpleCyclicGraph;

	ConstrainedBFSVisitorTest() {

		add_edge(0, 1, simpleAcyclicGraph);
		add_edge(0, 2, simpleAcyclicGraph);
		add_edge(2, 3, simpleAcyclicGraph);

		add_edge(0, 1, simpleCyclicGraph);
		add_edge(1, 3, simpleCyclicGraph);
		add_edge(0, 2, simpleCyclicGraph);
		add_edge(2, 3, simpleCyclicGraph);

	}
};

TEST_F(ConstrainedBFSVisitorTest, IdentifyUniquePath)
{
	int start = 0;
	int goal = 3;
	int minDepth = 0;
	int maxDepth = 2;
	int maxBranches = 3;

	DefaultColorMap<Graph> colorMap;
	ConstrainedBFSVisitor<Graph> visitor(start, goal, minDepth, maxDepth, maxBranches, colorMap);
	breadthFirstSearch(start, simpleAcyclicGraph, colorMap, visitor);

	AllPathsSearchResult<V> result = visitor.uniquePathToGoal();

	ASSERT_EQ(result.resultCode, FOUND_PATH);
	ASSERT_EQ(result.paths.size(), 1u);
	ASSERT_EQ(result.paths.at(0).str(), "0,2,3");
}

TEST_F(ConstrainedBFSVisitorTest, RespectMaxDepthLimit)
{
	int start = 0;
	int goal = 3;
	int minDepth = 0;
	int maxDepth = 1;
	int maxBranches = 3;

	DefaultColorMap<Graph> colorMap;
	ConstrainedBFSVisitor<Graph> visitor(start, goal, minDepth, maxDepth, maxBranches, colorMap);
	breadthFirstSearch(start, simpleAcyclicGraph, colorMap, visitor);

	EXPECT_EQ(visitor.getMaxDepthVisited(), 1);
	EXPECT_EQ(visitor.uniquePathToGoal().resultCode, NO_PATH);
}

TEST_F(ConstrainedBFSVisitorTest, RespectMinDepthLimit)
{
	int start = 0;
	int goal = 3;
	int minDepth = 3;
	int maxDepth = 10;
	int maxBranches = 3;

	DefaultColorMap<Graph> colorMap;
	ConstrainedBFSVisitor<Graph> visitor(start, goal, minDepth, maxDepth, maxBranches, colorMap);
	breadthFirstSearch(start, simpleAcyclicGraph, colorMap, visitor);

	EXPECT_EQ(visitor.uniquePathToGoal().resultCode, NO_PATH);
}

TEST_F(ConstrainedBFSVisitorTest, IdentifyMultiplePaths)
{
	int start = 0;
	int goal = 3;
	int minDepth = 0;
	int maxDepth = 3;
	int maxBranches = 3;

	DefaultColorMap<Graph> colorMap;
	ConstrainedBFSVisitor<Graph> visitor(start, goal, minDepth, maxDepth, maxBranches, colorMap);
	breadthFirstSearch(start, simpleCyclicGraph, colorMap, visitor);

	EXPECT_EQ(visitor.uniquePathToGoal().resultCode, TOO_MANY_PATHS);
}

TEST_F(ConstrainedBFSVisitorTest, ReturnMultiplePaths)
{
	int start = 0;
	int goal = 3;
	int minDepth = 0;
	int maxDepth = 3;
	int maxBranches = 3;

	DefaultColorMap<Graph> colorMap;
	ConstrainedBFSVisitor<Graph> visitor(start, goal, minDepth, maxDepth, maxBranches, colorMap);
	breadthFirstSearch(start, simpleCyclicGraph, colorMap, visitor);

	AllPathsSearchResult<V> result = visitor.pathsToGoal(2);

	EXPECT_EQ(result.resultCode, FOUND_PATH);
	ASSERT_EQ(result.paths.size(), 2u);

	string path1 = result.paths[0].str();
	string path2 = result.paths[1].str();

	EXPECT_TRUE(path1 != path2);
	ASSERT_TRUE(path1 == "0,1,3" || path1 == "0,2,3");
	ASSERT_TRUE(path2 == "0,1,3" || path2 == "0,2,3");
}

TEST_F(ConstrainedBFSVisitorTest, RespectBranchLimit)
{
	int start = 0;
	int goal = 3;
	int minDepth = 0;
	int maxDepth = 3;
	int maxBranches = 1;

	DefaultColorMap<Graph> colorMap;
	ConstrainedBFSVisitor<Graph> visitor(start, goal, minDepth, maxDepth, maxBranches, colorMap);
	breadthFirstSearch(start, simpleAcyclicGraph, colorMap, visitor);

	EXPECT_EQ(visitor.uniquePathToGoal().resultCode, TOO_MANY_BRANCHES);
}

}
