#ifndef OPTIONS_H
#define OPTIONS_H

#include <string>

namespace opt {
	/** k-mer length */
	extern unsigned kmerSize;

	/** read length */
	extern unsigned readLen;

	/** input FASTA path */
	extern std::string fastaFile;

	/** Parse the specified command line. */
	void parse(int argc, char* const* argv);
}

#endif
