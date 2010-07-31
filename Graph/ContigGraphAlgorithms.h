#ifndef CONTIGGRAPHALGORITHMS_H
#define CONTIGGRAPHALGORITHMS_H 1

/** Assemble an unambiguous path starting at vertex v. */
template<typename Graph, typename OutIt>
OutIt assemble(const Graph& g, typename Graph::vertex_descriptor v,
		OutIt out)
{
	typedef typename Graph::vertex_descriptor vertex_descriptor;
	for (; g.contiguous_out(v); v = g.vertex(g[v].front().target()))
		*out++ = v;
	*out++ = v;
	return out;
}

#endif
