#include "Graph/HashGraph.h"
#include <boost/tuple/tuple.hpp> // for boost::tie
#include <gtest/gtest.h>

using namespace std;

TEST(HashGraphTest, AddEdge)
{
	HashGraph<string> g;

	string a("a");
	string b("b");

	HashGraph<string>::edge_descriptor edge(a, b);
	
	add_edge(a, b, g);
	
	HashGraph<string>::out_edge_iterator ei, ei_end;
	boost::tie(ei, ei_end) = out_edges(a, g);

	ASSERT_EQ(edge, *ei);
}
