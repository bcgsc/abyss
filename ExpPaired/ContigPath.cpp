#include "ContigPath.h"
#include <algorithm>
#include <istream>
#include <iterator>
#include <ostream>

using namespace std;

void ContigPath::prependPath(const ContigPath& other)
{
	m_path.insert(m_path.begin(), other.m_path.begin(), other.m_path.end());
}

//
//
//
void ContigPath::appendPath(const ContigPath& other)
{
	m_path.insert(m_path.end(), other.m_path.begin(), other.m_path.end());
}

//
// Reverse the path and flip every merge node
//
void ContigPath::reverse(bool flipNodes)
{ 
	std::reverse(m_path.begin(), m_path.end());
	
	if(flipNodes)
	{
		size_t maxIdx = getNumNodes();
		for(size_t idx = 0; idx < maxIdx; ++idx)
		{
			getNode(idx).flip();
		}
	}
}

//
// Extract a subset
//
ContigPath ContigPath::extractNodes(size_t start, size_t end)
{
	ContigPath np;
	for(; start < end; ++start)
	{
		np.appendNode(getNode(start));
	}
	return np;
}

//
// Write the path to the stream
//
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
