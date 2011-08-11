#ifndef FASTA_INDEX_H
#define FASTA_INDEX_H 1

#include "IOUtil.h"
#include <algorithm>
#include <cassert>
#include <fstream>
#include <iterator> // for ostream_iterator
#include <limits> // for numeric_limits
#include <string>
#include <utility>
#include <vector>

/** A record of an indexed FASTA file. */
struct FAIRecord
{
	size_t offset;
	size_t size;
	std::string id;

	FAIRecord() : offset(0), size(0) { }
	FAIRecord(size_t offset, size_t size, const std::string& id)
		: offset(offset), size(size), id(id) { }

	friend std::ostream& operator<<(std::ostream& out,
			const FAIRecord& o)
	{
		return out << o.id << '\t' << o.size << '\t' << o.offset
			<< '\t' << o.size << '\t' << o.size + 1;
	}

	friend std::istream& operator>>(std::istream& in,
			FAIRecord& o)
	{
		size_t lineLen, lineBinLen;
		in >> o.id >> o.size >> o.offset >> lineLen >> lineBinLen;
		if (!in)
			return in;
		assert(o.size == lineLen || lineLen == lineBinLen);
		return in.ignore(std::numeric_limits<std::streamsize>::max(),
				'\n');
	}
};

/** An indexed FASTA (fai) file. */
class FastaIndex
{
	struct CompareOffset
	{
		/** Used with upper_bound. */
		bool operator()(size_t a, const FAIRecord& b) const
		{
			return a < b.offset;
		}
	};

  public:
	/** Index the specified FASTA file. */
	void index(const std::string& path)
	{
		m_data.clear();
		std::ifstream in(path.c_str());
		assert_good(in, path);
		char c;
		for (std::string id; in >> c && getline(in, id);) {
			assert(c == '>');
			size_t i = id.find(' ');
			if (i != std::string::npos)
				id.erase(i);
			assert(!id.empty());
			std::streampos offset = in.tellg();
			in.ignore(std::numeric_limits<std::streamsize>::max(),
					'\n');
			size_t n = in.gcount();
			assert(n > 0);
			m_data.push_back(FAIRecord(offset, n - 1, id));
		}
		assert(in.eof());
	}

	/** Translate a file offset to a sequence:position coordinate. */
	std::pair<std::string, size_t> operator[](size_t offset) const
	{
		Data::const_iterator it = std::upper_bound(
				m_data.begin(), m_data.end(),
				offset, CompareOffset());
		assert(it != m_data.begin());
		--it;
		assert(it != m_data.end());
		assert(it->offset <= offset);
		assert(offset < it->offset + it->size);
		return std::make_pair(it->id, offset - it->offset);
	}

	/** Write this index to a stream in SAM format. */
	void writeSAMHeader(std::ostream& out) const
	{
		for (Data::const_iterator it = m_data.begin();
				it != m_data.end(); ++it)
			out << "@SQ\tSN:" << it->id
				<< "\tLN:" << it->size << '\n';
	}

	/** Write this index to a stream. */
	friend std::ostream& operator<<(std::ostream& out,
			const FastaIndex& o)
	{
		std::copy(o.m_data.begin(), o.m_data.end(),
				std::ostream_iterator<FAIRecord>(out, "\n"));
		return out;
	}

	/** Read a FASTA index from a stream. */
	friend std::istream& operator>>(std::istream& in,
			FastaIndex& o)
	{
		assert(in.good());
		o.m_data.clear();
		for (FAIRecord rec; in >> rec;) {
			if (!o.m_data.empty())
				assert(rec.offset > o.m_data.back().offset);
			o.m_data.push_back(rec);
		}
		assert(in.eof());
		assert(!o.m_data.empty());
		return in;
	}

  private:
	typedef std::vector<FAIRecord> Data;
	Data m_data;
};

#endif
