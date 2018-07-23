#ifndef FASTAREADER_H
#define FASTAREADER_H 1

#include "Common/Sequence.h"
#include "Common/StringUtil.h" // for chomp
#include <cassert>
#include <cstdlib> // for exit
#include <fstream>
#include <istream>
#include <limits> // for numeric_limits
#include <ostream>
#include <cstdio>

/** Read a FASTA, FASTQ, export, qseq or SAM file. */
class FastaReader {
	public:
		enum {
			/** Fold lower-case characters to upper-case. */
			FOLD_CASE = 0, NO_FOLD_CASE = 1,
			/** Convert to standard quality. */
			NO_CONVERT_QUALITY = 0, CONVERT_QUALITY = 2,
		};
		bool flagFoldCase() { return ~m_flags & NO_FOLD_CASE; }
		bool flagConvertQual() { return m_flags & CONVERT_QUALITY; }

		FastaReader(const char* path, int flags, int len = 0);

		~FastaReader();

		Sequence read(std::string& id, std::string& comment,
				char& anchor, std::string& qual);

		/** Split the fasta file into nsections and seek to the start
		 * of section. */
		void split(unsigned section, unsigned nsections);

		/** Return whether this stream is at end-of-file. */
		bool eof() const {
			if (is_fasta)
				return m_file_state == -1;
			else
				return m_in->eof();
		};

		/** Return true if failbit or badbit of stream is set. */
		bool fail() const {
			if (is_fasta)
				return m_file_state < -1;
			else
				return m_in->fail();
		};

		/** Return whether this stream is good. */
		operator const void*() const {
			if (is_fasta) {
				return m_file_state >= 0 ? this : NULL;
			} else {
				return *m_in ? this : NULL;
			}
		}

		/** Return the next character of this stream. */
		int peek() {
			if (is_fasta)
				if (!eof())
					if (q)
						return '@';
					else
						return '>';
				else
					return EOF;
			else
				return m_in->peek();
		 }

		/** Interface for manipulators. */
		/*FastaReader& operator>>(std::istream& (*f)(std::istream&))
		{
			f(m_in);
			return *this;
		}*/

		/** Returns the number of unchaste reads. */
		unsigned unchaste() const { return m_unchaste; }

		FastaReader& operator >>(Sequence& seq)
		{
			std::string id, comment, qual;
			char anchor;
			seq = this->read(id, comment, anchor, qual);
			return *this;
		}

	private:
		/** Read a single line. */
		std::istream& getline(std::string& s)
		{
			if (std::getline(*m_in, s)) {
				chomp(s, '\r');
				m_line++;
			}
			return *m_in;
		}

		/** Ignore the specified number of lines. */
		std::istream& ignoreLines(unsigned n)
		{
			for (unsigned i = 0; i < n; ++i) {
				if (m_in->ignore(
						std::numeric_limits<std::streamsize>::max(),
						'\n'))
					m_line++;
			}
			return *m_in;
		}

		std::ostream& die();
		bool isChaste(const std::string& s, const std::string& line);
		void checkSeqQual(const std::string& s, const std::string& q);

		const char* m_path;
		std::ifstream* m_fin;
		std::istream* m_in;

		bool is_fasta = false; // Determines what method of reading to use
		bool q = false;
		void* m_file;
		int m_file_state = -2;
		void* seq_data;

		/** Flags indicating parsing options. */
		int m_flags;

		/** Number of lines read. */
		unsigned m_line;

		/** Count of unchaste reads. */
		unsigned m_unchaste;

		/** Position of the end of the current section. */
		std::streampos m_end;
		long m_file_end;

		/** Trim sequences to this length. 0 is unlimited. */
		const int m_maxLength;
};

/** A FASTA record. */
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

	FastaRecord() : anchor(0) { }
	FastaRecord(const std::string& id, const std::string& comment,
			const Sequence& seq)
		: id(id), comment(comment), anchor(0), seq(seq) { }

	operator Sequence() const { return seq; }

	FastaRecord& operator=(const std::string& s)
	{
		seq = s;
		return *this;
	}

	size_t size() const { return seq.size(); }

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
		if (o.qual.empty())
			return out << static_cast<const FastaRecord&>(o);
		out << '@' << o.id;
		if (!o.comment.empty())
			out << ' ' << o.comment;
		return out << '\n' << o.seq << "\n"
			"+\n" << o.qual << '\n';
	}
};

#endif //FASTAREADER_H
