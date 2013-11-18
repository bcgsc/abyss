#ifndef FASTAWRITER_H
#define FASTAWRITER_H 1

#include "Sequence.h"
#include <cstdio>

/** Output a FASTA file. */
class FastaWriter {
	public:
		// Constructor opens file
		FastaWriter(const char* path, bool append = false);

		// Destructor closes it
		~FastaWriter();

		/** Write a sequence with a comment. */
		void WriteSequence(const Sequence& seq, unsigned id,
				unsigned multiplicity, const std::string& comment);

		/** Write a sequence. */
		void WriteSequence(const Sequence& seq, unsigned id,
				unsigned multiplicity)
		{
			WriteSequence(seq, id, multiplicity, "");
		}

		void WriteSequence(const Sequence& seq, unsigned long long id,
				const std::string& comment);

		void WriteSequence(const Sequence& seq, const std::string& id,
				const std::string& comment);

	private:
		const char *m_path;
		FILE* m_fileHandle;
};

#endif
