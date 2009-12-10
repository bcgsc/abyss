#ifndef CONTIGPATH_H
#define CONTIGPATH_H 1

#include "Dictionary.h"
#include "PairUtils.h" // for LinearNumKey
#include <cassert>
#include <istream>
#include <ostream>
#include <sstream>
#include <vector>

struct MergeNode
{
	LinearNumKey id;
	bool isRC;
	
	void flip() { isRC = (isRC) ? 0 : 1; }

	bool operator ==(const MergeNode& o) const
	{
		return id == o.id && isRC == o.isRC;
	}

	bool operator <(const MergeNode& o) const
	{
		return id != o.id ? id < o.id : isRC < o.isRC;
	}

	friend std::ostream& operator<<(std::ostream& out,
			const MergeNode& o)
	{
		return out << g_contigIDs.key(o.id) << ',' << o.isRC;
	}

	friend std::istream& operator>>(std::istream& in, MergeNode& o)
	{
		char c;
		unsigned cid;
		in >> cid >> c >> o.isRC;
		std::stringstream s;
		s << cid;
		o.id = g_contigIDs.serial(s.str());
		if (in.good())
			assert(c == ',');
		return in;
	}
};

class ContigPath : public std::vector<MergeNode>
{
	public:
		ContigPath::ContigPath() { }

		void reverse(bool flipNodes);

		/** Return a subsequence of this path. */
		ContigPath extractNodes(size_t start, size_t end)
		{
			std::vector<MergeNode> v(&(*this)[start], &(*this)[end]);
			return ContigPath(v);
		}

		friend std::ostream& operator<<(std::ostream& out,
				const ContigPath& object);
		friend std::istream& operator>>(std::istream& in,
				ContigPath& object);

	private:
		explicit ContigPath(const std::vector<MergeNode>& v)
			: std::vector<MergeNode>(v) { }
};

#endif
