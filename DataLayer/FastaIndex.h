#ifndef FASTA_INDEX_H
#define FASTA_INDEX_H 1

#include "IOUtil.h"
#include <algorithm>
#include <cassert>
#include <fstream>
#include <limits> // for numeric_limits
#include <string>
#include <utility>
#include <vector>

/** An indexed FASTA (fai) file. */
class FastaIndex
{
	typedef std::vector<std::pair<size_t, std::string> > Index;

	struct CompareFirst
	{
		/** Used with upper_bound. */
		bool operator()(size_t a, const Index::value_type& b) const
		{
			return a < b.first;
		}
	};

  public:
	FastaIndex(const std::string& path)
	{
		std::ifstream in(path.c_str());
		assert_good(in, path);
		std::string id;
		size_t len, pos, lineLen, lineBinLen;
		while (in >> id >> len >> pos >> lineLen >> lineBinLen) {
			in.ignore(std::numeric_limits<std::streamsize>::max(),
					'\n');
			assert(len == lineLen || lineLen == lineBinLen);
			if (!m_data.empty())
				assert(pos > m_data.back().first);
			m_data.push_back(Index::value_type(pos, id));
		}
		assert(in.eof());
		assert(!m_data.empty());
	}

	/** Translate a file offset to a sequence:position coordinate. */
	std::pair<std::string, size_t> operator[](size_t pos) const
	{
		Index::const_iterator it = std::upper_bound(
				m_data.begin(), m_data.end(),
				pos, CompareFirst());
		assert(it != m_data.begin());
		--it;
		assert(it != m_data.end());
		assert(it->first <= pos);
		return std::make_pair(it->second, pos - it->first);
	}

  private:
	Index m_data;
};

#endif
