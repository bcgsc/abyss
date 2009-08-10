#ifndef FASTAWRITER_H
#define FASTAWRITER_H 1

#include "Sequence.h"
#include <cstdio>

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

	private:
		const char *m_path;
		FILE* m_fileHandle;
};

#endif
