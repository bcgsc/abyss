#include "ContigPath.h"
#include <algorithm>
#include <istream>
#include <ostream>
#include <iterator>

//
//
//
ContigPath::ContigPath()
{
	
}

//
// Append a single node to the list
//
void ContigPath::appendNode(const MergeNode& mn)
{
	m_path.push_back(mn);
}

//
//
//
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
//
//
size_t ContigPath::findFirstOf(LinearNumKey id)
{
	size_t maxIdx = getNumNodes();
	for(size_t idx = 0; idx < maxIdx; ++idx)
	{
		if(getNode(idx).id == id)
		{
			return idx;
		}
	}
	
	return maxIdx;
}

//
//
//
size_t ContigPath::findLastOf(LinearNumKey id)
{
	size_t maxIdx = getNumNodes();
	for(size_t idx = 0; idx < maxIdx; ++idx)
	{
		size_t tIdx = maxIdx - idx - 1;
		if(getNode(tIdx).id == id)
		{
			return tIdx;
		}
	}
	
	return maxIdx;
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
std::ostream& operator<<(std::ostream& out, const ContigPath& object)
{
	std::copy(object.m_path.begin(), object.m_path.end(), std::ostream_iterator<MergeNode>(out, " "));
	return out;
}

//
// Read the path from the stream
//
std::istream& operator>>(std::istream& in, ContigPath& object)
{
	std::copy(std::istream_iterator<MergeNode>(in), std::istream_iterator<MergeNode>(), back_inserter(object.m_path));
	return in;
}
