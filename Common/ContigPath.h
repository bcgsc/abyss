#ifndef CONTIGPATH_H
#define CONTIGPATH_H 1

#include "ContigNode.h"
#include <algorithm>
#include <cassert>
#include <functional>
#include <istream>
#include <sstream>
#include <string>
#include <ostream>
#include <vector>

/** A sequence of ContigNode. */
typedef std::vector<ContigNode> ContigPath;

/** Reverse and complement the specified path. */
template<typename T>
void reverseComplement(T first, T last)
{
	std::reverse(first, last);
	for (T it = first; it < last; ++it)
		if (!it->ambiguous())
			*it ^= 1;
}

/** Return the reverse complement of the specified path. */
static inline ContigPath reverseComplement(const ContigPath& path)
{
	ContigPath rc(path.rbegin(), path.rend());
	for (ContigPath::iterator it = rc.begin(); it < rc.end(); ++it)
		if (!it->ambiguous())
			*it ^= 1;
	return rc;
}

static inline std::ostream& operator<<(std::ostream& out,
		const ContigPath& o)
{
	assert(!o.empty());
	ContigPath::const_iterator it = o.begin();
	out << get(g_contigNames, *it);
	for (++it; it != o.end(); ++it)
		out << ' ' << get(g_contigNames, *it);
	return out;
}

static inline std::istream& operator>>(std::istream& in,
		ContigPath& o)
{
	o.clear();
	std::string s;
	if (getline(in, s)) {
		std::istringstream ss(s);
		for (std::string name; ss >> name;)
			o.push_back(find_vertex(name, g_contigNames));
		assert(ss.eof());
	}
	return in;
}

static inline void markSeenInPath(std::istream& in,
		std::vector<bool>& marked)
{
	assert(in.good());
	std::string s;
	ContigPath path;
	while (in >> s >> path) {
		if (path.empty()) {
			size_t i = get(g_contigNames, s);
			assert(i < marked.size());
			marked[i] = true;
		}
		for (ContigPath::const_iterator it = path.begin();
				it != path.end(); ++it) {
			if (!it->ambiguous()) {
				size_t i = it->contigIndex();
				assert(i < marked.size());
				marked[i] = true;
			}
		}
	}
	assert(in.eof());
}

#endif
