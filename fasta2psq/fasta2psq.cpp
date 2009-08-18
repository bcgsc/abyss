#include "fasta2psqOptions.h"
#include "FastaReader.h"
#include "PackedSeqWriter.h"
#include <cstdio>

int main(int argc, char* const* argv)
{
	pp_opt::parse(argc, argv);

	printf("Reading from file: %s\n", pp_opt::fastaFile.c_str());
	printf("Writing k-mer to %s\n", pp_opt::outFile.c_str());

	FastaReader reader(pp_opt::fastaFile.c_str());
	PackedSeqWriter writer(pp_opt::outFile.c_str());

	for (Sequence seq; reader >> seq;)
		writer.WriteSequence(seq);
	assert(reader.eof());

	return 0;
}
