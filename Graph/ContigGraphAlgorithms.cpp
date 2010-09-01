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

typedef ContigGraph<DirectedGraph<ContigProperties> > Graph;

/** Assemble unambiguous paths. Write the paths to out. */
void assemble(Graph& g, ostream& out)
{
	typedef Graph::vertex_descriptor vertex_descriptor;
	typedef Graph::vertex_iterator vertex_iterator;
	pair<vertex_iterator, vertex_iterator> vit = g.vertices();
	for (vertex_iterator v = vit.first; v != vit.second; ++v) {
		if (!contiguous_out(g, *v) || contiguous_in(g, *v)
				|| is_palindrome(g, *out_edges(*v, g).first))
			continue;
		vector<vertex_descriptor> path;
		assemble(g, *v, back_inserter(path));
		assert(path.size() >= 2);
		assert(path.front() != path.back());
		out << ContigID::create() << '\t' << ContigPath(path) << '\n';
		merge(g, path.begin(), path.end());
		remove_vertex_if(g, path.begin(), path.end(),
				not1(std::mem_fun_ref(&ContigNode::ambiguous)));
	}
}
