#ifndef CONTIGPATH_H
#define CONTIGPATH_H 1

#include "ContigNode.h"
#include "StringUtil.h"
#include <algorithm>
#include <cassert>
#include <functional>
#include <istream>
#include <iterator>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>

typedef ContigNode MergeNode;

class ContigPath : public std::vector<MergeNode>
{
	typedef std::vector<MergeNode> Vector;

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
					std::mem_fun_ref(&MergeNode::flip));
		}
};

std::ostream& operator<<(std::ostream& out, const ContigPath& o)
{
	assert(!o.empty());
	ContigPath::const_iterator last = o.end() - 1;
	copy(o.begin(), last, std::ostream_iterator<MergeNode>(out, " "));
	return out << *last;
}

std::istream& operator>>(std::istream& in, ContigPath& o)
{
	copy(std::istream_iterator<MergeNode>(in),
			std::istream_iterator<MergeNode>(),
			back_inserter(o));
	return in;
}

#endif
