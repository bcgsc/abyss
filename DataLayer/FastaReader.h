#ifndef FASTAREADER_H
#define FASTAREADER_H 1

#include "Sequence.h"
#include <cassert>
#include <fstream>
#include <istream>
#include <ostream>

class FastaReader {
	public:
		enum {
			/** Discard any sequence containing an N. */
			DISCARD_N = 0, KEEP_N = 1,
			/** Fold lower-case characters to upper-case. */
			FOLD_CASE = 0, NO_FOLD_CASE = 2,
			/** Convert to standard quality. */
			NO_CONVERT_QUALITY = 0, CONVERT_QUALITY = 4,
		};
		bool flagDiscardN() { return ~m_flags & KEEP_N; }
		bool flagFoldCase() { return ~m_flags & NO_FOLD_CASE; }
		bool flagConvertQual() { return m_flags & CONVERT_QUALITY; }

		FastaReader(const char* path, int flags);
		~FastaReader();

		Sequence read(std::string& id, std::string& comment,
				char& anchor, std::string& qual);

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
			std::string id, comment, qual;
			char anchor;
			seq = this->read(id, comment, anchor, qual);
			return *this;
		}

	private:
		const char* m_inPath;
		std::ifstream m_inFile;
		std::istream& m_fileHandle;

		/** Flags indicating parsing options. */
		int m_flags;

		/** Count of unchaste reads. */
		unsigned m_unchaste;

		/** Count of non-ACGT reads. */
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

	FastaRecord() { }
	FastaRecord(const std::string& id, const std::string& comment,
			const Sequence& seq)
		: id(id), comment(comment), anchor(0), seq(seq) { }

	friend FastaReader& operator >>(FastaReader& in, FastaRecord& o)
	{
		std::string q;
		o.seq = in.read(o.id, o.comment, o.anchor, q);
		return in;
	}

	friend std::ostream& operator <<(std::ostream& out,
			const FastaRecord& o)
	{
		out << '>' << o.id;
		if (!o.comment.empty())
			out << ' ' << o.comment;
		return out << '\n' << o.seq << '\n';
	}
};

/** A FASTQ record. */
struct FastqRecord : FastaRecord
{
	/** Quality */
	std::string qual;

	FastqRecord() { }
	FastqRecord(const std::string& id, const std::string& comment,
			const Sequence& seq, const std::string& qual)
		: FastaRecord(id, comment, seq), qual(qual)
	{
		assert(seq.length() == qual.length());
	}

	friend FastaReader& operator >>(FastaReader& in, FastqRecord& o)
	{
		o.seq = in.read(o.id, o.comment, o.anchor, o.qual);
		return in;
	}

	friend std::ostream& operator <<(std::ostream& out,
			const FastqRecord& o)
	{
		out << '@' << o.id;
		if (!o.comment.empty())
			out << ' ' << o.comment;
		return out << '\n' << o.seq << "\n"
			"+\n" << o.qual << '\n';
	}
};

#endif //FASTAREADER_H
