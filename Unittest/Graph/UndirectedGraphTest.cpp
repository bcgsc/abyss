#include "Common/Hash.h"
#include "Graph/UndirectedGraph.h"
#include "Graph/Properties.h"
#include <boost/tuple/tuple.hpp>
#include <gtest/gtest.h>
#include <string>
#include <sstream>

using namespace std;

typedef UndirectedGraph<NoProperty, NoProperty> Graph;
typedef boost::graph_traits<Graph>::vertex_descriptor V;
typedef boost::graph_traits<Graph>::edge_descriptor E;
typedef boost::graph_traits<Graph>::edge_iterator EdgeIterator;

TEST(UndirectedGraphTest, edgeComparison)
{
	Graph g;
	V u = add_vertex(NoProperty(), g);
	V v = add_vertex(NoProperty(), g);
	E uv = E(u, v);
	E vu = E(v, u);
	ASSERT_EQ(uv, vu);
	ASSERT_EQ(uv.hashCode(), vu.hashCode());
}

TEST(UndirectedGraphTest, edges)
{
	Graph g;

	V u = add_vertex(NoProperty(), g);
	V v = add_vertex(NoProperty(), g);

	/*
	 * Internally, both edge orientations (u, v) and (v, u) are
	 * added to the graph data structures.  Make sure that
	 * this isn't exposed externally and that we only see one edge.
	 */

	E e;
	bool inserted;
	boost::tie(e, inserted) = add_edge(u, v, NoProperty(), g);
	ASSERT_TRUE(inserted);
	ASSERT_EQ(1u, num_edges(g));

	/* make sure iterator only visits one version of the edge */

	EdgeIterator eit, eit_end;
	boost::tie(eit, eit_end) = edges(g);
	ASSERT_NE(eit, eit_end);
	ASSERT_EQ(e, *eit);
	++eit;
	ASSERT_EQ(eit, eit_end);
}
