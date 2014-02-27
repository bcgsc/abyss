#ifndef BIDIRECTION_BFS_VISITOR_H
#define BIDIRECTION_BFS_VISITOR_H 1

#include "Graph/Path.h"
#include "Common/Warnings.h"

enum BFSVisitorResult { SUCCESS, ABORT_SEARCH, SKIP_ELEMENT };

template <class Graph>
class BidirectionalBFSVisitor {

public:

	typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
	typedef typename boost::graph_traits<Graph>::edge_descriptor Edge;

	BidirectionalBFSVisitor() { }
	virtual ~BidirectionalBFSVisitor() { }

	virtual BFSVisitorResult discover_vertex(const Vertex& u, const Graph& g, Direction dir, unsigned numFrontierNodes)
	{
		SUPPRESS_UNUSED_WARNING(u);
		SUPPRESS_UNUSED_WARNING(g);
		SUPPRESS_UNUSED_WARNING(dir);
		SUPPRESS_UNUSED_WARNING(numFrontierNodes);
		return SUCCESS;
	}

	virtual void examine_vertex(const Vertex& u, const Graph& g, Direction dir)
	{
		SUPPRESS_UNUSED_WARNING(u);
		SUPPRESS_UNUSED_WARNING(g);
		SUPPRESS_UNUSED_WARNING(dir);
	}

	virtual void finish_vertex(const Vertex& u, const Graph& g, Direction dir)
	{
		SUPPRESS_UNUSED_WARNING(u);
		SUPPRESS_UNUSED_WARNING(g);
		SUPPRESS_UNUSED_WARNING(dir);
	}

	virtual void examine_edge(const Edge& e, const Graph& g, Direction dir)
	{
		SUPPRESS_UNUSED_WARNING(e);
		SUPPRESS_UNUSED_WARNING(g);
		SUPPRESS_UNUSED_WARNING(dir);
	}

	virtual BFSVisitorResult common_edge(const Edge& e, const Graph& g, Direction dir)
	{
		SUPPRESS_UNUSED_WARNING(e);
		SUPPRESS_UNUSED_WARNING(g);
		SUPPRESS_UNUSED_WARNING(dir);
		return SUCCESS;
	}

	virtual BFSVisitorResult tree_edge(const Edge& e, const Graph& g, Direction dir)
	{
		SUPPRESS_UNUSED_WARNING(e);
		SUPPRESS_UNUSED_WARNING(g);
		SUPPRESS_UNUSED_WARNING(dir);
		return SUCCESS;
	}

	virtual BFSVisitorResult non_tree_edge(const Edge& e, const Graph& g, Direction dir)
	{
		SUPPRESS_UNUSED_WARNING(e);
		SUPPRESS_UNUSED_WARNING(g);
		SUPPRESS_UNUSED_WARNING(dir);
		return SUCCESS;
	}

	virtual void gray_target(const Edge& e, const Graph& g, Direction dir)
	{
		SUPPRESS_UNUSED_WARNING(e);
		SUPPRESS_UNUSED_WARNING(g);
		SUPPRESS_UNUSED_WARNING(dir);
	}

	virtual void black_target(const Edge& e, const Graph& g, Direction dir)
	{
		SUPPRESS_UNUSED_WARNING(e);
		SUPPRESS_UNUSED_WARNING(g);
		SUPPRESS_UNUSED_WARNING(dir);
	}

};

#endif

