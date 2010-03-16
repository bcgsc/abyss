#ifndef CONTIGPATH_H
#define CONTIGPATH_H 1

#include "ContigNode.h"
#include <algorithm>
#include <cassert>
#include <functional>
#include <istream>
#include <iterator>
#include <ostream>
#include <vector>

class ContigPath : public std::vector<ContigNode>
{
	typedef std::vector<ContigNode> Vector;

	public:
		ContigPath() { }

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

		/** The separator to print between ContigNode. */
		static const char* separator;
};
const char* ContigPath::separator;

std::ostream& operator<<(std::ostream& out, const ContigPath& o)
{
	const char* separator = o.separator == NULL ? " " : o.separator;
	assert(!o.empty());
	ContigPath::const_iterator last = o.end() - 1;
	copy(o.begin(), last,
			std::ostream_iterator<ContigNode>(out, separator));
	return out << *last;
}

std::istream& operator>>(std::istream& in, ContigPath& o)
{
	copy(std::istream_iterator<ContigNode>(in),
			std::istream_iterator<ContigNode>(),
			back_inserter(o));
	return in;
}

#endif
