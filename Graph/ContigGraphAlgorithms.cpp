/**
 * Contig graph algorithms.
 * @author Shaun Jackman <sjackman@bcgsc.ca>.
 */

#include "ContigGraphAlgorithms.h"
#include "ContigPath.h"
#include <cassert>
#include <functional>
#include <vector>

using namespace std;

typedef ContigGraph<ContigProperties> Graph;
typedef Graph::vertex_descriptor vertex_descriptor;
typedef Graph::vertex_iterator vertex_iterator;
typedef vector<vertex_descriptor> Path;

/** Return the contig properties of the specified path. */
static ContigProperties calculateProperties(const Graph& g,
		const Path& path)
{
	ContigProperties vp(opt::k - 1, 0);
	for (Path::const_iterator it = path.begin();
			it != path.end(); ++it) {
		assert(!it->ambiguous());
		vp += g[*it];
	}
	return vp;
}

/** Merge the specified path and update the graph g. */
static void mergePath(Graph& g, const Path& path, ostream& out)
{
	vertex_descriptor v = g.add_vertex(calculateProperties(g, path));
	copy_in_edges(g, path.front(), v);
	copy_out_edges(g, path.back(), v);
	for_each(path.begin(), path.end(),
			bind1st(mem_fun(&Graph::clear_vertex), &g));
	for_each(path.begin(), path.end(),
			bind1st(mem_fun(&Graph::remove_vertex), &g));

	ContigID id = ContigID::create();
	assert(ContigID(v) == id);
	out << id << '\t' << ContigPath(path) << '\n';
}

/** Assemble unambiguous paths. Write the paths to out. */
void assemble(Graph& g, ostream& out)
{
	pair<vertex_iterator, vertex_iterator> vit = g.vertices();
	for (vertex_iterator v = vit.first; v != vit.second; ++v) {
		if (contiguous_in(g, *v) || !contiguous_out(g, *v))
			continue;
		Path path;
		assemble(g, *v, back_inserter(path));
		assert(path.size() >= 3);
		assert(path.front() != path.back());
		// Output only the canonical path.
		if (path.front() < path.back())
			mergePath(g, path, out);
	}
}
