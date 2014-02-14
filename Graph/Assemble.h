#ifndef ASSEMBLE_H
#define ASSEMBLE_H 1

#include "Common/ContigNode.h" // for ContigIndexMap
#include "Common/Iterator.h" // for output_iterator_traits
#include "Graph/DepthFirstSearch.h"

using boost::graph_traits;

/** Return true if this edge is contiguous. */
template <typename Graph>
bool isContiguous(const Graph &g,
		typename graph_traits<Graph>::edge_descriptor e)
{
	return out_degree(source(e, g), g) == 1
		&& in_degree(target(e, g), g) == 1;
}

/** Assemble contigous paths. Write the paths to out. */
template <typename OutIt>
class AssembleVisitor : public boost::default_dfs_visitor
{
  public:
	AssembleVisitor(OutIt it) : m_it(it) { }

	template <typename Vertex, typename Graph>
	void discover_vertex(const Vertex& u, Graph&)
	{
		m_path.push_back(u);
	}

	template <typename Vertex, typename Graph>
	void finish_vertex(const Vertex&, Graph&)
	{
		finishContig();
	}

	template <typename Edge, typename Graph>
	void examine_edge(const Edge& e, const Graph& g)
	{
		if (!isContiguous(g, e))
			finishContig();
	}

  private:
	void finishContig()
	{
		if (m_path.size() > 1)
			*m_it++ = m_path;
		m_path.clear();
	}

	OutIt m_it;
	typename output_iterator_traits<OutIt>::value_type m_path;
};

/** Assemble unambiguous paths. Write the paths to out. */
template<typename Graph, typename OutIt>
void assembleDFS(const Graph& g, OutIt out, bool ss = false)
{
	(void)ss;
	typedef boost::vector_property_map<
		boost::default_color_type, ContigIndexMap> ColorMap;
	depthFirstSearch(g, AssembleVisitor<OutIt>(out),
			ColorMap(num_vertices(g) / 2), ss);
}

#endif
