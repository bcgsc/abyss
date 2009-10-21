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

class ContigPath
{
	public:
		ContigPath() { }

		/** Append a single node to this path. */
		void appendNode(const MergeNode& mn) { m_path.push_back(mn); }

		/** Prepend the specified path to this path. */
		void prependPath(const ContigPath& o)
		{
			m_path.insert(m_path.begin(),
					o.m_path.begin(), o.m_path.end());
		}

		/** Append the specified path to this path. */
		void appendPath(const ContigPath& o)
		{
			m_path.insert(m_path.end(),
					o.m_path.begin(), o.m_path.end());
		}

		// Get the number of nodes in the path
		size_t getNumNodes() const { return m_path.size(); }
		
		// Get the node with a specified index
		MergeNode& getNode(size_t idx)
		{
			assert(idx < m_path.size());
			return m_path[idx];
		}
		
		// Get the node with a specified index
		const MergeNode& getNode(size_t idx) const
		{
			assert(idx < m_path.size());
			return m_path[idx];
		}

		// reverse the path
		void reverse(bool flipNodes);

		/** Return a subsequence of this path. */
		ContigPath extractNodes(size_t start, size_t end)
		{
			std::vector<MergeNode> v(&m_path[start], &m_path[end]);
			return ContigPath(v);
		}

		bool operator <(const ContigPath& o) const
		{
			return m_path < o.m_path;
		}

		// Insertion/Extraction operators
		friend std::ostream& operator<<(std::ostream& out, const ContigPath& object);
		friend std::istream& operator>>(std::istream& in, ContigPath& object);

	private:
		explicit ContigPath(const std::vector<MergeNode>& v)
			: m_path(v) { }

		std::vector<MergeNode> m_path;
};

#endif
