#ifndef FASTACONCAT_H
#define FASTACONCAT_H 1

#include "FastaReader.h"
#include <cassert>
#include <vector>
#include <iostream>

class FastaConcat 
{

  public:

	FastaConcat(char** first, char** last, int flags)
		: m_fail(false), m_flags(flags), m_fileIndex(0),
		m_reader(NULL)
	{
		assert(first != last);
		for (char** p = first; p != last; p++)
			m_filenames.push_back(*p);
		m_reader = new FastaReader(
			m_filenames.at(m_fileIndex).c_str(), flags);
	}

	~FastaConcat()
	{
		if(m_reader != NULL)
			delete m_reader;
	}

	bool eof() const
	{
		return m_fileIndex == m_filenames.size();
	}

	/** Return true if all of the streams are still good. */
	operator void*() const
	{
		return m_fail ? NULL : const_cast<FastaConcat*>(this);
	}

	template <typename Record>
	friend FastaConcat& operator>>(
			FastaConcat& in, Record& o)
	{
		in.m_fail = false;
		for (; in.m_fileIndex < in.m_filenames.size(); in.m_fileIndex++) {
			assert(in.m_reader != NULL);
			*in.m_reader >> o;
			if (!in.m_reader->fail()) {
				return in;
			} else if (in.m_reader->fail() && !in.m_reader->eof()) {
				in.m_fail = true;
				return in;
			} else {
				delete in.m_reader;
				in.m_reader = NULL;
				if (in.m_fileIndex < in.m_filenames.size() - 1) {
					in.m_reader = new FastaReader(
						in.m_filenames.at(in.m_fileIndex + 1).c_str(),
						in.m_flags);
				}
			}
		}
		// set fail when attempting to read at eof
		// (like normal iostreams)
		in.m_fail = true;
		return in;
	}

  private:

	/** Emulates failbit of iostream */
	bool m_fail;

	/** FastaReader flags */
	int m_flags;

	/** Index of current file */
	unsigned m_fileIndex;

	/** FastaReader for current file */
	FastaReader* m_reader;

	/** List of filenames to read */
	std::vector<std::string> m_filenames;
};

#endif
