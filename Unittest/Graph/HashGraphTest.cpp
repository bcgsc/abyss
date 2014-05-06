#include "Graph/HashGraph.h"
#include <string>
#include <gtest/gtest.h>

using namespace std;

namespace {

class HashGraphTest : public ::testing::Test {

protected:

	typedef HashGraph<string> Graph;
	typedef boost::graph_traits<Graph>::edge_descriptor edge_descriptor;
	typedef boost::graph_traits<Graph>::out_edge_iterator out_edge_iterator;
	typedef boost::graph_traits<Graph>::in_edge_iterator in_edge_iterator;
	typedef boost::graph_traits<Graph>::vertex_descriptor vertex_descriptor;
	typedef boost::graph_traits<Graph>::vertex_iterator vertex_iterator;

	Graph simpleCyclicGraph;
	
	string a;
	string b;
	string c;
	string d;

	unordered_set<string> vertexSet;

	HashGraphTest() : a("a"), b("b"), c("c"), d("d") {

		add_edge(a, b, simpleCyclicGraph);
		add_edge(a, c, simpleCyclicGraph);
		add_edge(b, d, simpleCyclicGraph);
		add_edge(c, d, simpleCyclicGraph);

		vertexSet.insert(a);
		vertexSet.insert(b);
		vertexSet.insert(c);
		vertexSet.insert(d);

	}
};

TEST_F(HashGraphTest, out_edge_iterator)
{
	edge_descriptor expectedEdge1(a, b);
	edge_descriptor expectedEdge2(a, c);
	
	out_edge_iterator ei, ei_end;
	boost::tie(ei, ei_end) = out_edges(a, simpleCyclicGraph);

	ASSERT_TRUE(ei != ei_end);
	edge_descriptor edge1 = *ei;
	EXPECT_TRUE(edge1 == expectedEdge1 || edge1 == expectedEdge2);

	ei++;
	ASSERT_TRUE(ei != ei_end);
	edge_descriptor edge2 = *ei;
	EXPECT_TRUE(edge2 != edge1);
	EXPECT_TRUE(edge2 == expectedEdge1 || edge2 == expectedEdge2);
	
	ei++;
	EXPECT_TRUE(ei == ei_end);
}

TEST_F(HashGraphTest, vertex_iterator)
{
	vertex_iterator vi, vi_begin, vi_end;
	boost::tie(vi_begin, vi_end) = vertices(simpleCyclicGraph);
	ASSERT_TRUE(vi_begin != vi_end);

	unsigned count = 0;
	vertex_descriptor v;
	vi = vi_begin;
	for(; vi != vi_end; vi++, count++) {
		// check that current vertex is different from
		// previous vertex
		if (vi != vi_begin)
			EXPECT_TRUE(*vi != v);
		vertex_descriptor v = *vi;
		EXPECT_TRUE(vertexSet.find(v) != vertexSet.end());
	}

	EXPECT_EQ(4U, count);
}

}
