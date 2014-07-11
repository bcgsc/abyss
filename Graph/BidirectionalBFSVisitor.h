#ifndef BIDIRECTION_BFS_VISITOR_H
#define BIDIRECTION_BFS_VISITOR_H 1

#include "Graph/Path.h"

enum BFSVisitorResult { SUCCESS, ABORT_SEARCH, SKIP_ELEMENT };

template <class Graph>
class BidirectionalBFSVisitor {

public:

	typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
	typedef typename boost::graph_traits<Graph>::edge_descriptor Edge;
	typedef unsigned NumFrontierNodes;

	BidirectionalBFSVisitor() { }
	virtual ~BidirectionalBFSVisitor() { }

	virtual BFSVisitorResult discover_vertex(const Vertex&, const Graph&, Direction, NumFrontierNodes)
	{
		return SUCCESS;
	}

	virtual void examine_vertex(const Vertex&, const Graph&, Direction)
	{
	}

	virtual void finish_vertex(const Vertex&, const Graph&, Direction)
	{
	}

	virtual void examine_edge(const Edge&, const Graph&, Direction)
	{
	}

	virtual BFSVisitorResult common_edge(const Edge&, const Graph&, Direction)
	{
		return SUCCESS;
	}

	virtual BFSVisitorResult tree_edge(const Edge&, const Graph&, Direction)
	{
		return SUCCESS;
	}

	virtual BFSVisitorResult non_tree_edge(const Edge&, const Graph&, Direction)
	{
		return SUCCESS;
	}

	virtual void gray_target(const Edge&, const Graph&, Direction)
	{
	}

	virtual void black_target(const Edge&, const Graph&, Direction)
	{
	}

};

#endif

