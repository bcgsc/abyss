#include "Graph/Path.h"
#include "Graph/ConstrainedBidiBFSVisitor.h"
#include "Graph/BidirectionalBFS.h"
#include <boost/graph/adjacency_list.hpp>
#include <boost/property_map/property_map.hpp>
#include <limits>
#include <gtest/gtest.h>

using namespace std;
using namespace boost;

typedef adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS> Graph;
typedef graph_traits<Graph>::vertex_descriptor V;
typedef graph_traits<Graph>::edge_descriptor e;
typedef std::vector< Path<V> > PathList;

static const size_t NO_MEM_LIMIT = std::numeric_limits<std::size_t>::max();

namespace {

class ConstrainedBidiBFSVisitorTest : public ::testing::Test {

protected:

	Graph simpleAcyclicGraph;
	Graph simpleCyclicGraph;
	Graph cyclicGraph;

	ConstrainedBidiBFSVisitorTest() {

		add_edge(0, 1, simpleAcyclicGraph);
		add_edge(0, 2, simpleAcyclicGraph);
		add_edge(2, 3, simpleAcyclicGraph);

		add_edge(0, 1, simpleCyclicGraph);
		add_edge(1, 3, simpleCyclicGraph);
		add_edge(0, 2, simpleCyclicGraph);
		add_edge(2, 3, simpleCyclicGraph);

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

TEST_F(ConstrainedBidiBFSVisitorTest, IdentifyUniquePath)
{
	ConstrainedBidiBFSVisitor<Graph>
		visitor(simpleAcyclicGraph, 0, 3, 1, 1, 3, 2, NO_MEM_LIMIT);
	bidirectionalBFS(simpleAcyclicGraph, 0, 3, visitor);

	Path<V> uniquePath;
	PathSearchResult result = visitor.uniquePathToGoal(uniquePath);

	ASSERT_EQ(FOUND_PATH, result);
	ASSERT_EQ("0,2,3", uniquePath.str());
}

TEST_F(ConstrainedBidiBFSVisitorTest, StartEqualsGoal)
{
	ConstrainedBidiBFSVisitor<Graph>
		visitor(simpleAcyclicGraph, 0, 0, 1, 1, 1, 2, NO_MEM_LIMIT);
	bidirectionalBFS(simpleAcyclicGraph, 0, 0, visitor);

	Path<V> uniquePath;
	PathSearchResult result = visitor.uniquePathToGoal(uniquePath);

	ASSERT_EQ(FOUND_PATH, result);
	ASSERT_EQ("0", uniquePath.str());
}

TEST_F(ConstrainedBidiBFSVisitorTest, SingleEdgeToGoal)
{
	ConstrainedBidiBFSVisitor<Graph>
		visitor(simpleAcyclicGraph, 0, 1, 1, 1, 2, 2, NO_MEM_LIMIT);
	bidirectionalBFS(simpleAcyclicGraph, 0, 1, visitor);

	Path<V> uniquePath;
	PathSearchResult result = visitor.uniquePathToGoal(uniquePath);

	ASSERT_EQ(FOUND_PATH, result);
	ASSERT_EQ("0,1", uniquePath.str());
}

TEST_F(ConstrainedBidiBFSVisitorTest, RespectMaxPathLength)
{
	ConstrainedBidiBFSVisitor<Graph>
		visitor(cyclicGraph, 0, 6, 4, 5, 6, 2, NO_MEM_LIMIT);
	bidirectionalBFS(cyclicGraph, 0, 6, visitor);

	vector< Path<V> > paths;
	PathSearchResult result = visitor.pathsToGoal(paths);

	// We expect the fourth path ("0,1,2,3,4,5,6")
	// to be excluded by the max path length limit (6).

	set<string> expectedPaths;
	expectedPaths.insert("0,1,3,5,6");
	expectedPaths.insert("0,1,2,3,5,6");
	expectedPaths.insert("0,1,3,4,5,6");

	EXPECT_EQ(FOUND_PATH, result);
	ASSERT_EQ(3u, paths.size());

	// check that each path is one of the expected ones
	for (unsigned i = 0; i < 3; i++)
		ASSERT_TRUE(expectedPaths.find(paths[i].str()) != expectedPaths.end());

	// check that each path is unique
	for (unsigned i = 0; i < 2; i++)
		ASSERT_TRUE(paths[i].str() != paths[i+1].str());
}

TEST_F(ConstrainedBidiBFSVisitorTest, RespectMinPathLength)
{
	ConstrainedBidiBFSVisitor<Graph>
		visitor(cyclicGraph, 0, 6, 4, 6, 7, 2, NO_MEM_LIMIT);
	bidirectionalBFS(cyclicGraph, 0, 6, visitor);

	vector< Path<V> > paths;
	PathSearchResult result = visitor.pathsToGoal(paths);

	// We expect the first path ("0,1,3,5,6")
	// to be excluded by the min path length limit (6).

	set<string> expectedPaths;
	expectedPaths.insert("0,1,2,3,5,6");
	expectedPaths.insert("0,1,3,4,5,6");
	expectedPaths.insert("0,1,2,3,4,5,6");

	EXPECT_EQ(FOUND_PATH, result);
	ASSERT_EQ(3u, paths.size());

	// check that each path is one of the expected ones
	for (unsigned i = 0; i < 3; i++)
		ASSERT_TRUE(expectedPaths.find(paths[i].str()) != expectedPaths.end());

	// check that each path is unique
	for (unsigned i = 0; i < 2; i++)
		ASSERT_TRUE(paths[i].str() != paths[i+1].str());
}

TEST_F(ConstrainedBidiBFSVisitorTest, RespectMaxPathsLimit)
{
	ConstrainedBidiBFSVisitor<Graph>
		visitor(simpleCyclicGraph, 0, 3, 1, 1, 3, 2, NO_MEM_LIMIT);
	bidirectionalBFS(simpleCyclicGraph, 0, 3, visitor);

	Path<V> uniquePath;
	EXPECT_EQ(visitor.uniquePathToGoal(uniquePath), TOO_MANY_PATHS);
}

TEST_F(ConstrainedBidiBFSVisitorTest, ReturnMultiplePaths)
{
	ConstrainedBidiBFSVisitor<Graph>
		visitor(simpleCyclicGraph, 0, 3, 2, 1, 3, 2, NO_MEM_LIMIT);
	bidirectionalBFS(simpleCyclicGraph, 0, 3, visitor);

	PathList paths;
	PathSearchResult result = visitor.pathsToGoal(paths);

	EXPECT_EQ(FOUND_PATH, result);
	ASSERT_EQ(2u, paths.size());

	string path1 = paths[0].str();
	string path2 = paths[1].str();

	EXPECT_TRUE(path1 != path2);
	ASSERT_TRUE(path1 == "0,1,3" || path1 == "0,2,3");
	ASSERT_TRUE(path2 == "0,1,3" || path2 == "0,2,3");
}

TEST_F(ConstrainedBidiBFSVisitorTest, RespectMaxBranches)
{
	ConstrainedBidiBFSVisitor<Graph>
		visitor(simpleCyclicGraph, 0, 3, 2, 1, 3, 1, NO_MEM_LIMIT);
	bidirectionalBFS(simpleCyclicGraph, 0, 3, visitor);

	PathList paths;
	PathSearchResult result = visitor.pathsToGoal(paths);

	EXPECT_EQ(TOO_MANY_BRANCHES, result);
	EXPECT_EQ(0u, paths.size());
	// expect early exit from traversal
	EXPECT_EQ(3u, visitor.getNumNodesVisited());
}

TEST_F(ConstrainedBidiBFSVisitorTest, NoLimitForBranches)
{
	ConstrainedBidiBFSVisitor<Graph>
		visitor(simpleCyclicGraph, 0, 3, 2, 1, 3, NO_LIMIT, NO_MEM_LIMIT);
	bidirectionalBFS(simpleCyclicGraph, 0, 3, visitor);

	PathList paths;
	PathSearchResult result = visitor.pathsToGoal(paths);

	EXPECT_EQ(FOUND_PATH, result);
	ASSERT_EQ(2u, paths.size());

	string path1 = paths[0].str();
	string path2 = paths[1].str();

	EXPECT_TRUE(path1 != path2);
	ASSERT_TRUE(path1 == "0,1,3" || path1 == "0,2,3");
	ASSERT_TRUE(path2 == "0,1,3" || path2 == "0,2,3");
}

}
