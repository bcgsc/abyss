#include "Graph/Path.h"
#include "Graph/ExtendPath.h"
#include <boost/graph/adjacency_list.hpp>
#include <gtest/gtest.h>

using namespace std;

typedef boost::adjacency_list<boost::vecS, boost::vecS,
	boost::bidirectionalS> Graph;

// note: vertex_descriptor for adjacency_list<> is int
typedef boost::graph_traits<Graph>::vertex_descriptor Vertex;
// note: edge_descriptor for adjacency_list<> is int
typedef boost::graph_traits<Graph>::edge_descriptor Edge;

TEST(extendPath, lookAhead)
{
	unsigned depth;

	/* case 1: simple path */

	/*
	 * 0--1--2
	 */

	Graph g1;
	add_edge(0, 1, g1);
	add_edge(1, 2, g1);

	depth = 1;
	ASSERT_TRUE(lookAhead(1, FORWARD, depth, g1));
	ASSERT_TRUE(lookAhead(1, REVERSE, depth, g1));
	ASSERT_FALSE(lookAhead(2, FORWARD, depth, g1));
	ASSERT_FALSE(lookAhead(0, REVERSE, depth, g1));

	depth = 2;
	ASSERT_FALSE(lookAhead(1, FORWARD, depth, g1));
	ASSERT_FALSE(lookAhead(1, REVERSE, depth, g1));
	ASSERT_TRUE(lookAhead(0, FORWARD, depth, g1));
	ASSERT_TRUE(lookAhead(2, REVERSE, depth, g1));

	/* case 2: with branching */

	/*
	 *      2
	 *     /
	 * 0--1
	 *     \
	 *      3--4
	 */

	Graph g2;
	add_edge(0, 1, g2);
	add_edge(1, 2, g2);
	add_edge(1, 3, g2);
	add_edge(3, 4, g2);

	depth = 3;
	ASSERT_TRUE(lookAhead(0, FORWARD, depth, g2));

	depth = 4;
	ASSERT_FALSE(lookAhead(0, FORWARD, depth, g2));
}

TEST(extendPath, noExtension)
{
	// Graph containing a single edge.

	Graph g;
	add_edge(0, 1, g);
	Path<Vertex> path;
	path.push_back(0);
	path.push_back(1);

	extendPath(path, FORWARD, g);
	ASSERT_EQ(2, path.size());

	extendPath(path, REVERSE, g);
	ASSERT_EQ(2, path.size());
}

TEST(extendPath, extendForward)
{
	/*
	 *      2
	 *     /
	 * 0--1
	 *     \
	 *      3
	 */

	Graph g;
	add_edge(0, 1, g);
	add_edge(1, 2, g);
	add_edge(1, 3, g);

	Path<Vertex> expectedPath;
	expectedPath.push_back(0);
	expectedPath.push_back(1);

	Path<Vertex> path;
	path.push_back(0);
	ASSERT_EQ(1, path.size());

	extendPath(path, FORWARD, g);
	ASSERT_EQ(2, path.size());
	ASSERT_EQ(expectedPath, path);
}

TEST(extendPath, extendReverse)
{
	/*
	 *  0
	 *   \
	 *    2--3
	 *   /
	 *  1
	 */

	Graph g;
	add_edge(0, 2, g);
	add_edge(1, 2, g);
	add_edge(2, 3, g);

	Path<Vertex> expectedPath;
	expectedPath.push_back(2);
	expectedPath.push_back(3);

	Path<Vertex> path;
	path.push_back(3);
	ASSERT_EQ(1, path.size());

	extendPath(path, REVERSE, g);
	ASSERT_EQ(2, path.size());
	ASSERT_EQ(expectedPath, path);
}

TEST(extendPath, bidirectional)
{
	/*
	 *  0         5
	 *   \       /
	 *    2--3--4
	 *   /       \
	 *  1         6
	 */

	Graph g;
	add_edge(0, 2, g);
	add_edge(1, 2, g);
	add_edge(2, 3, g);
	add_edge(3, 4, g);
	add_edge(4, 5, g);
	add_edge(4, 6, g);

	Path<Vertex> expectedPath;
	expectedPath.push_back(2);
	expectedPath.push_back(3);
	expectedPath.push_back(4);

	Path<Vertex> path;
	path.push_back(3);
	ASSERT_EQ(1, path.size());

	extendPath(path, FORWARD, g);
	extendPath(path, REVERSE, g);
	EXPECT_EQ(3, path.size());
	ASSERT_EQ(expectedPath, path);
}

TEST(extendPath, withTrimming)
{
	const unsigned trimLen = 1;

	/*
	 *       2
	 *      /
	 *  0--1--3--4
	 */

	Graph g;
	add_edge(0, 1, g);
	add_edge(1, 2, g);
	add_edge(1, 3, g);
	add_edge(3, 4, g);

	Path<Vertex> expectedPath;
	expectedPath.push_back(0);
	expectedPath.push_back(1);
	expectedPath.push_back(3);
	expectedPath.push_back(4);

	Path<Vertex> path;
	path.push_back(0);

	extendPath(path, FORWARD, g, trimLen);
	ASSERT_EQ(4, path.size());
	ASSERT_EQ(expectedPath, path);

	/*
	 *       2  4
	 *      /  /
	 *  0--1--3
	 *         \
	 *          5
	 */

	Graph g2;
	add_edge(0, 1, g2);
	add_edge(1, 2, g2);
	add_edge(1, 3, g2);
	add_edge(3, 4, g2);
	add_edge(3, 5, g2);

	Path<Vertex> expectedPath2;
	expectedPath2.push_back(0);
	expectedPath2.push_back(1);
	expectedPath2.push_back(3);

	Path<Vertex> path2;
	path2.push_back(0);

	extendPath(path2, FORWARD, g2, trimLen);
	EXPECT_EQ(3, path2.size());
	ASSERT_EQ(expectedPath2, path2);
}
