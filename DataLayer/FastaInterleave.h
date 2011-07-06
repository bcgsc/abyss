#ifndef FASTAINTERLEAVE_H
#define FASTAINTERLEAVE_H 1

#include "FastaReader.h"
#include <cassert>
#include <vector>

class FastaInterleave {
	typedef FastaReader Stream;
	typedef std::vector<Stream*> Streams;

  public:
	FastaInterleave(char** first, char** last, int flags)
		: m_it(m_streams.begin()), m_fail(false)
	{
		assert(first != last);
		m_streams.reserve(last - first);
		for (char** p = first; p < last; ++p)
			m_streams.push_back(new Stream(*p, flags));
		m_it = m_streams.begin();
	}

	~FastaInterleave()
	{
		for (Streams::iterator it = m_streams.begin();
				it != m_streams.end(); ++it)
			delete *it;
	}

	/** Return true if all the streams are eof. */
	bool eof() const
	{
		for (Streams::const_iterator it = m_streams.begin();
				it != m_streams.end(); ++it)
			if (!(*it)->eof())
				return false;
		return true;
	}

	/** Return true if any of the streams are good. */
	operator void*() const
	{
		return m_fail ? NULL : const_cast<FastaInterleave*>(this);
	}

	/** Extract one record from the next stream. */
	template <typename Record>
	friend FastaInterleave& operator>>(
			FastaInterleave& in, Record& o)
	{
		for (unsigned i = 0; i < in.m_streams.size(); ++i) {
			assert(in.m_it != in.m_streams.end());
			bool good = **in.m_it >> o;
			if (++in.m_it == in.m_streams.end())
				in.m_it = in.m_streams.begin();
			if (good) {
				in.m_fail = false;
				return in;
			}
		}
		in.m_fail = true;
		return in;
	}

  private:
	/** The streams. */
	Streams m_streams;

	/** The next stream from which to read. */
	Streams::iterator m_it;

	/** True when all streams have failed. */
	bool m_fail;
};

#endif
