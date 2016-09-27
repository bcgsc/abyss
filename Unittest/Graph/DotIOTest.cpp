#include "Common/ContigID.h"
#include "Common/Dictionary.h"
#include "Graph/GraphIO.h"
#include "Graph/DotIO.h"
#include "Graph/UndirectedGraph.h"
#include "Graph/DirectedGraph.h"
#include "Graph/Properties.h"
#include <boost/tuple/tuple.hpp>
#include <gtest/gtest.h>
#include <string>
#include <sstream>

using namespace std;

typedef UndirectedGraph<NoProperty, NoProperty> UGraph;
typedef boost::graph_traits<UGraph>::vertex_descriptor UV;
typedef boost::graph_traits<UGraph>::edge_descriptor UE;
typedef boost::graph_traits<UGraph>::edge_iterator UEIt;

typedef DirectedGraph<NoProperty, NoProperty> DGraph;
typedef boost::graph_traits<DGraph>::vertex_descriptor DV;
typedef boost::graph_traits<DGraph>::edge_descriptor DE;
typedef boost::graph_traits<DGraph>::edge_iterator DEIt;

TEST(DotIOTest, UndirectedGraph)
{
	/* reset state of global contig ID map */
	g_contigNames.clear();

	string graphviz;
	stringstream ss;

	/* undirected graph, unquoted IDs */

	graphviz =
		"graph {\n"
		"     a;\n"
		"     b;\n"
		"     a -- b;\n"
		"}\n";

	ss << graphviz;
	ASSERT_TRUE(ss.good());

	UGraph g;
	read_dot(ss, g, DisallowParallelEdges());
	ASSERT_TRUE(!ss.fail());

	/* check that correct graph structure was created */

	ASSERT_EQ(2u, num_vertices(g));
	UV a = find_vertex("a", g);
	ASSERT_LT(a, num_vertices(g));
	UV b = find_vertex("b", g);
	ASSERT_LT(b, num_vertices(g));

	ASSERT_EQ(1u, num_edges(g));
	UEIt ueit, ueit_end;
	boost::tie(ueit, ueit_end) = edges(g);
	ASSERT_NE(ueit, ueit_end);
	UE e = UE(a, b);
	ASSERT_EQ(e, *ueit);
	ASSERT_TRUE(source(*ueit, g) == a || target(*ueit, g) == a);
	ASSERT_TRUE(source(*ueit, g) == b || target(*ueit, g) == b);
	++ueit;
	ASSERT_EQ(ueit, ueit_end);

}

TEST(DotIOTest, DirectedGraph)
{
	/* reset state of global contig ID map */
	g_contigNames.clear();

	string graphviz;
	stringstream ss;

	/* directed graph, quoted IDs */

	graphviz =
		"digraph {\n"
		"     \"1+\";\n"
		"     \"1-\";\n"
		"     \"2+\";\n"
		"     \"2-\";\n"
		"     \"1+\" -> \"2-\";\n"
		"}\n";

	ss << graphviz;
	ASSERT_TRUE(ss.good());

	DGraph g;
	read_dot(ss, g, DisallowParallelEdges());
	ASSERT_TRUE(!ss.fail());

	/* check that correct graph structure was created */

	ASSERT_EQ(4u, num_vertices(g));
	DV v1 = find_vertex("1+", g);
	ASSERT_LT(v1.index(), num_vertices(g));
	DV v2 = find_vertex("2-", g);
	ASSERT_LT(v2.index(), num_vertices(g));

	ASSERT_EQ(1u, num_edges(g));
	DEIt ueit, ueit_end;
	boost::tie(ueit, ueit_end) = edges(g);
	ASSERT_NE(ueit, ueit_end);
	DE e = DE(v1, v2);
	ASSERT_EQ(e, *ueit);
	ASSERT_TRUE(source(*ueit, g) == v1);
	ASSERT_TRUE(target(*ueit, g) == v2);
	++ueit;
	ASSERT_EQ(ueit, ueit_end);

}
