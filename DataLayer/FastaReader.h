#ifndef FASTAREADER_H
#define FASTAREADER_H

#include "Sequence.h"
#include <fstream>
#include <istream>

class FastaReader {
	public:
		FastaReader(const char* path, bool discardN = true);
		enum { KEEP_N = false, DISCARD_N = true };

		~FastaReader();

		Sequence read(std::string& id, std::string& comment,
				char& anchor);

		/** Return whether this stream is at end-of-file. */
		bool eof() const { return m_fileHandle.eof(); };

		/** Return whether this stream is good. */
		operator void*() const { return m_fileHandle; }

		/** Returns the number of unchaste reads. */
		unsigned unchaste() const { return m_unchaste; }

		/** Returns the number of reads containing non-ACGT
		 * characters. */
		unsigned nonACGT() const { return m_nonacgt; }

		FastaReader& operator >>(Sequence& seq)
		{
			std::string id;
			std::string comment;
			char anchor;
			seq = this->read(id, comment, anchor);
			return *this;
		}

	private:
		const char* m_inPath;
		std::ifstream m_inFile;
		std::istream& m_fileHandle;
		bool m_discardN;
		unsigned m_unchaste;
		unsigned m_nonacgt;
};

struct FastaRecord
{
	/** Identifier */
	std::string id;
	/** Comment following the first white-space of the header */
	std::string comment;
	/** Anchor base for a colour-space sequence */
	char anchor;
	/** The sequence */
	Sequence seq;

	friend FastaReader& operator >>(FastaReader& in, FastaRecord& o)
	{
		o.seq = in.read(o.id, o.comment, o.anchor);
		return in;
	}
};

#endif //FASTAREADER_H
