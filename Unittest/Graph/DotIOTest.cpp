#include "Graph/GraphIO.h"
#include "Graph/DotIO.h"
#include "Graph/UndirectedGraph.h"
#include "Graph/Properties.h"
#include <gtest/gtest.h>
#include <string>
#include <sstream>

using namespace std;

TEST(DotIOTest, read_dot)
{
	string undirectedGraph =
		"graph {\n"
		"     \"a\";\n"
		"     \"b\";\n"
		"     \"a\" -- \"b\";\n"
		"}\n";

	stringstream ss;
	ss << undirectedGraph;
	ASSERT_TRUE(ss.good());

	UndirectedGraph<NoProperty, NoProperty> g;
	read_dot(ss, g, DisallowParallelEdges());
	ASSERT_TRUE(!ss.fail());
}
