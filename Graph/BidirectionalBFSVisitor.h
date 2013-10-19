#ifndef BIDIRECTION_BFS_VISITOR_H
#define BIDIRECTION_BFS_VISITOR_H 1

#include "Common/Warnings.h"

template <class Graph>
class BidirectionalBFSVisitor {

public:

	typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
	typedef typename boost::graph_traits<Graph>::edge_descriptor Edge;

	BidirectionalBFSVisitor() { }

	virtual void discover_vertex(Vertex u, const Graph& g, Direction dir)
	{
		SUPPRESS_UNUSED_WARNING(u);
		SUPPRESS_UNUSED_WARNING(g);
		SUPPRESS_UNUSED_WARNING(dir);
	}

	virtual void examine_vertex(Vertex u, const Graph& g, Direction dir)
	{
		SUPPRESS_UNUSED_WARNING(u);
		SUPPRESS_UNUSED_WARNING(g);
		SUPPRESS_UNUSED_WARNING(dir);
	}

	virtual void finish_vertex(Vertex u, const Graph& g, Direction dir)
	{
		SUPPRESS_UNUSED_WARNING(u);
		SUPPRESS_UNUSED_WARNING(g);
		SUPPRESS_UNUSED_WARNING(dir);
	}

	virtual void examine_edge(Edge e, const Graph& g, Direction dir)
	{
		SUPPRESS_UNUSED_WARNING(e);
		SUPPRESS_UNUSED_WARNING(g);
		SUPPRESS_UNUSED_WARNING(dir);
	}

	virtual void common_edge(Edge e, const Graph& g)
	{
		SUPPRESS_UNUSED_WARNING(e);
		SUPPRESS_UNUSED_WARNING(g);
	}

	virtual void tree_edge(Edge e, const Graph& g, Direction dir)
	{
		SUPPRESS_UNUSED_WARNING(e);
		SUPPRESS_UNUSED_WARNING(g);
		SUPPRESS_UNUSED_WARNING(dir);
	}

	virtual void non_tree_edge(Edge e, const Graph& g, Direction dir)
	{
		SUPPRESS_UNUSED_WARNING(e);
		SUPPRESS_UNUSED_WARNING(g);
		SUPPRESS_UNUSED_WARNING(dir);
	}

	virtual void gray_target(Edge e, const Graph& g, Direction dir)
	{
		SUPPRESS_UNUSED_WARNING(e);
		SUPPRESS_UNUSED_WARNING(g);
		SUPPRESS_UNUSED_WARNING(dir);
	}

	virtual void black_target(Edge e, const Graph& g, Direction dir)
	{
		SUPPRESS_UNUSED_WARNING(e);
		SUPPRESS_UNUSED_WARNING(g);
		SUPPRESS_UNUSED_WARNING(dir);
	}

};

#endif

