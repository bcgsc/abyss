#ifndef CONTIGPATH_H
#define CONTIGPATH_H 1

#include "ContigNode.h"
#include "ContigID.h"
#include <algorithm>
#include <cassert>
#include <functional>
#include <istream>
#include <iterator>
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
	std::for_each(first, last,
			std::mem_fun_ref(&ContigNode::flip));
}

/** Return the reverse complement of the specified path. */
static inline ContigPath reverseComplement(const ContigPath& path)
{
	ContigPath rc(path.rbegin(), path.rend());
	std::for_each(rc.begin(), rc.end(),
			std::mem_fun_ref(&ContigNode::flip));
	return rc;
}

static inline std::ostream& operator<<(std::ostream& out,
		const ContigPath& o)
{
	assert(!o.empty());
	ContigPath::const_iterator last = o.end() - 1;
	copy(o.begin(), last,
			std::ostream_iterator<ContigNode>(out, " "));
	return out << *last;
}

static inline std::istream& operator>>(std::istream& in,
		ContigPath& o)
{
	o.clear();
	std::string s;
	if (getline(in, s)) {
		std::istringstream ss(s);
		for (ContigNode u; ss >> u;)
			o.push_back(u);
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
			ContigID id(s);
			assert(marked.size() > id);
			marked[id] = true;
		}
		for (ContigPath::const_iterator it = path.begin();
				it != path.end(); ++it) {
			if (!it->ambiguous()) {
				ContigID id(*it);
				assert(marked.size() > id);
				marked[id] = true;
			}
		}
	}
	assert(in.eof());
}

#endif
