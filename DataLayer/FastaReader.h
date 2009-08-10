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

		// Read sequences into the vector as packed seqs
		// Returns true unless eof has been reached
		bool ReadSequences(SequenceVector& outseqs);

		// Returns true unless eof has been reached
		bool isGood()
		{
			return !m_fileHandle.eof() && m_fileHandle.peek() != EOF;
		}

		/** Return whether this stream is at end-of-file. */
		bool eof() const { return m_fileHandle.eof(); };

		/** Return whether this stream is good. */
		operator void*() const { return m_fileHandle; }

		// Returns the number of sequences containing non-ACGT
		// characters.
		unsigned getNonACGT() { return m_nonacgt; }

	private:
		const char* m_inPath;
		std::ifstream m_inFile;
		std::istream& m_fileHandle;
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
