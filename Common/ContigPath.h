#ifndef CONTIGPATH_H
#define CONTIGPATH_H 1

#include "ContigNode.h"
#include <algorithm>
#include <cassert>
#include <functional>
#include <istream>
#include <iterator>
#include <sstream>
#include <string>
#include <ostream>
#include <vector>

class ContigPath : public std::vector<ContigNode>
{
	typedef std::vector<ContigNode> Vector;

	public:
		ContigPath() { }
		explicit ContigPath(size_t n) : Vector(n) { }
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

		using std::vector<ContigNode>::erase;

		/** The separator to print between ContigNode. */
		static const char* separator;
};

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
			std::ostream_iterator<ContigNode>(out, o.separator));
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

#endif
