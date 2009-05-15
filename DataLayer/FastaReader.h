#ifndef FASTAREADER_H
#define FASTAREADER_H

#include "IFileReader.h"
#include "Sequence.h"
#include <fstream>
#include <istream>

class FastaReader : public IFileReader
{
	public:

		// Constructor opens file
		FastaReader(const char* path);

		// Destructor closes it
		~FastaReader();

		// Read a single sequence from the file
		Sequence ReadSequence(std::string& id);
		Sequence ReadSequence()
		{
			std::string id;
			return ReadSequence(id);
		}

		// Read sequences into the vector as packed seqs
		// Returns true unless eof has been reached
		virtual bool ReadSequences(SequenceVector& outseqs);
		
		// Returns true unless eof has been reached
		bool isGood();
				
		// Returns the number of sequences containing non-ACGT
		// characters.
		virtual unsigned getNonACGT() { return m_nonacgt; }

	private:
		const char* m_inPath;
		std::ifstream m_inFile;
		std::istream& m_fileHandle;
		unsigned m_nonacgt;
};

#endif //FASTAREADER_H
