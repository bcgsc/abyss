#include "Graph/Path.h"
#include "Graph/BreadthFirstSearch.h"
#include "Graph/ConstrainedBFSVisitor.h"
#include "Graph/DefaultColorMap.h"
#include <boost/graph/adjacency_list.hpp>
#include <boost/property_map/property_map.hpp>
#include <gtest/gtest.h>

using namespace std;
using namespace boost;

typedef adjacency_list<> Graph;
typedef ConstrainedBFSVisitor<Graph>::Path Path;
typedef ConstrainedBFSVisitor<Graph>::PathList PathList;

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
	breadthFirstSearch(simpleAcyclicGraph, start, visitor, colorMap);

	Path uniquePath;
	PathSearchResult result = visitor.uniquePathToGoal(uniquePath);

	ASSERT_EQ(result, FOUND_PATH);
	ASSERT_EQ(visitor.pathToString(uniquePath), "0,2,3");
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
	breadthFirstSearch(simpleAcyclicGraph, start, visitor, colorMap);

	Path uniquePath;

	EXPECT_EQ(visitor.getMaxDepthVisited(), 1);
	EXPECT_EQ(visitor.uniquePathToGoal(uniquePath), NO_PATH);
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
	breadthFirstSearch(simpleAcyclicGraph, start, visitor, colorMap);

	Path uniquePath;

	EXPECT_EQ(visitor.uniquePathToGoal(uniquePath), NO_PATH);
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
	breadthFirstSearch(simpleCyclicGraph, start, visitor, colorMap);

	Path uniquePath;

	EXPECT_EQ(visitor.uniquePathToGoal(uniquePath), TOO_MANY_PATHS);
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
	breadthFirstSearch(simpleCyclicGraph, start, visitor, colorMap);

	PathList paths;
	PathSearchResult result = visitor.pathsToGoal(paths, 2);

	EXPECT_EQ(result, FOUND_PATH);
	ASSERT_EQ(paths.size(), 2);

	string path1 = visitor.pathToString(paths[0]);
	string path2 = visitor.pathToString(paths[1]);

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
	breadthFirstSearch(simpleAcyclicGraph, start, visitor, colorMap);

	Path uniquePath;

	EXPECT_EQ(visitor.uniquePathToGoal(uniquePath), TOO_MANY_BRANCHES);
}

}

int main(int argc, char** argv)
{
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}






