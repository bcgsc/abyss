#ifndef FASTA_INDEX_H
#define FASTA_INDEX_H 1

#include "IOUtil.h"
#include <boost/tuple/tuple.hpp>
#include <algorithm>
#include <cassert>
#include <fstream>
#include <iterator> // for ostream_iterator
#include <string>
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
		return in >> Ignore('\n');
	}
};

/** An indexed FASTA (fai) file. */
class FastaIndex
{
	typedef std::vector<FAIRecord> Data;

	struct CompareOffset
	{
		/** Used with upper_bound. */
		bool operator()(size_t a, const FAIRecord& b) const
		{
			return a < b.offset;
		}
	};

  public:
	typedef boost::tuple<const FAIRecord&, size_t> SeqPos;

	/** Return the number of contigs. */
	size_t size() { return m_data.size(); }

	/** Return the size of the FASTA file. */
	size_t fileSize() const
	{
		assert(!m_data.empty());
		return m_data.back().offset + m_data.back().size + 1;
	}

	typedef Data::const_iterator const_iterator;
	const_iterator begin() const { return m_data.begin(); }
	const_iterator end() const { return m_data.end(); }

	/** Index the specified FASTA file. */
	void index(const std::string& path)
	{
		m_data.clear();
		std::ifstream in(path.c_str());
		assert_good(in, path);
		char c;
		for (std::string id;
				in >> c && in >> id && in >> Ignore('\n');) {
			assert(c == '>');
			assert(!id.empty());
			std::streampos offset = in.tellg();
			assert(offset > 0);
			in >> Ignore('\n');
			size_t n = in.gcount();
			assert(n > 0);
			m_data.push_back(FAIRecord(offset, n - 1, id));
		}
		assert(in.eof());
	}

	/** Translate a file offset to a sequence:position coordinate. */
	SeqPos operator[](size_t offset) const
	{
		Data::const_iterator it = std::upper_bound(
				m_data.begin(), m_data.end(),
				offset, CompareOffset());
		assert(it != m_data.begin());
		--it;
		assert(it != m_data.end());
		assert(it->offset <= offset);
		assert(offset < it->offset + it->size);
		return SeqPos(*it, offset - it->offset);
	}

	/** Write FASTA headers to the specified seekable stream. */
	void writeFASTAHeaders(std::ostream& out) const
	{
		assert(out);
		if (!out.seekp(0))
			return;
		for (Data::const_iterator it = m_data.begin();
				it != m_data.end(); ++it) {
			out << '>' << it->id << ' ';
			assert(it->offset > 0);
			if (!out.seekp(it->offset - 1))
				break;
			out << '\n';
			if (!out.seekp(it->offset + it->size))
				break;
			out << '\n';
		}
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

		// Count the number of records.
		if (in.seekg(0, std::ios::beg)) {
			size_t n = 0;
			while (in)
				if (in.get() == '\n')
					n++;
			o.m_data.reserve(n);
			in.clear();
			in.seekg(0, std::ios::beg);
			assert(in.good());
		} else
			in.clear();

		// Read the records.
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
	Data m_data;
};

#endif
