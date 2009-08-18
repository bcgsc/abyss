#ifndef FASTAREADER_H
#define FASTAREADER_H

#include "Sequence.h"
#include <fstream>
#include <istream>

class FastaReader {
	public:

		// Constructor opens file
		FastaReader(const char* path);

		// Destructor closes it
		~FastaReader();

		// Read a single sequence from the file
		Sequence ReadSequence(std::string& id, std::string& comment,
				char& anchor);
		Sequence ReadSequence(std::string& id, char& anchor)
		{
			std::string comment;
			return ReadSequence(id, comment, anchor);
		}
		Sequence ReadSequence(std::string& id)
		{
			char anchor;
			return ReadSequence(id, anchor);
		}
		Sequence ReadSequence()
		{
			std::string id;
			return ReadSequence(id);
		}

		// Returns true unless eof has been reached
		bool isGood()
		{
			return !m_fileHandle.eof() && m_fileHandle.peek() != EOF;
		}

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
			seq = this->ReadSequence();
			return *this;
		}

	private:
		const char* m_inPath;
		std::ifstream m_inFile;
		std::istream& m_fileHandle;
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
		o.seq = in.ReadSequence(o.id, o.comment, o.anchor);
		return in;
	}
};

#endif //FASTAREADER_H
