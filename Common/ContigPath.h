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
class ContigPath : public std::vector<ContigNode>
{
	typedef ContigNode T;
	typedef std::vector<T> Vector;

	public:
		ContigPath() { }
		explicit ContigPath(size_t n, const T& x = T())
			: Vector(n, x) { }
		explicit ContigPath(const Vector& v) : Vector(v) { }

		template <class InputIterator>
		ContigPath(InputIterator first, InputIterator last)
			: Vector(first, last) { }

		/** Reverse the path and flip every node. */
		void reverseComplement()
		{
			std::reverse(begin(), end());
			std::for_each(begin(), end(),
					std::mem_fun_ref(&ContigNode::flip));
		}

		using Vector::erase;
};

/** Return the reverse complement of the specified path. */
static inline ContigPath reverseComplement(const ContigPath& path)
{
	ContigPath rc(path.rbegin(), path.rend());
	std::for_each(rc.begin(), rc.end(),
			std::mem_fun_ref(&ContigNode::flip));
	return rc;
}

namespace std {
	template<>
	inline void swap(ContigPath& a, ContigPath& b) { a.swap(b); }
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
		copy(std::istream_iterator<ContigNode>(ss),
				std::istream_iterator<ContigNode>(),
				back_inserter(o));
		assert(ss.eof());
	}
	return in;
}

static inline void markSeenInPath(std::istream& in,
		std::vector<bool>& marked)
{
	assert(in.good());
	std::string id;
	ContigPath path;
	while (in >> id >> path) {
		assert(marked.size() > ContigID(id));
		if (path.empty())
			marked[ContigID(id)] = true;
		for (ContigPath::const_iterator it = path.begin();
				it != path.end(); ++it)
			if (!it->ambiguous())
				marked[ContigID(*it)] = true;
	}
	assert(in.eof());
}

#endif
