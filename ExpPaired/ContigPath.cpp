#include "ContigPath.h"
#include <algorithm>
#include <istream>
#include <iterator>
#include <ostream>

using namespace std;

/** Reverse the path and flip every node. */
void ContigPath::reverse(bool flipNodes)
{
	std::reverse(m_path.begin(), m_path.end());
	if (flipNodes) {
		size_t maxIdx = getNumNodes();
		for (size_t idx = 0; idx < maxIdx; ++idx)
			getNode(idx).flip();
	}
}

/** Write a path. */
ostream& operator<<(ostream& out, const ContigPath& object)
{
	vector<MergeNode>::const_iterator last = object.m_path.end() - 1;
	copy(object.m_path.begin(), last,
			ostream_iterator<MergeNode>(out, " "));
	return out << *last;
}

/** Read a path. */
istream& operator>>(istream& in, ContigPath& object)
{
	copy(istream_iterator<MergeNode>(in),
			istream_iterator<MergeNode>(),
			back_inserter(object.m_path));
	return in;
}
