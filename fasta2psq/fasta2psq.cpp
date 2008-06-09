#include "fasta2psq.h"
#include "FastaReader.h"
#include "PackedSeqWriter.h"
#include <cstdio>

int main(int argc, char* const* argv)
{	
	pp_opt::parse(argc, argv);
	
	printf("Writing sequences to %s \n", pp_opt::outFile.c_str());
	printf("Reading from file: %s\n", pp_opt::fastaFile.c_str());

	
	// write data
	FastaReader* reader = new FastaReader(pp_opt::fastaFile.c_str());
	PackedSeqWriter* writer = new PackedSeqWriter(pp_opt::outFile.c_str());
	
	bool stop = false;
	while(!stop)
	{
		SequenceVector seqs;
		stop = !reader->ReadSequences(seqs);
		for(SequenceVectorIterator iter = seqs.begin(); iter != seqs.end(); iter++)
		{
			writer->WriteSequence(*iter);
		}
	}
	
	delete reader;
	delete writer;

	return 0;
}
