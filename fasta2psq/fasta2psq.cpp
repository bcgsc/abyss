#include <stdio.h>

#include <vector>
#include <stdio.h>
#include <deque>
#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include "fasta2psq.h"
#include "FastaWriter.h"
#include "FastaReader.h"
#include "PackedSeqWriter.h"
#include "PackedSeqReader.h"


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
		PSequenceVector seqs;
		stop = !reader->ReadSequences(seqs);
		for(PSequenceVectorIterator iter = seqs.begin(); iter != seqs.end(); iter++)
		{
			writer->WriteSequence(*iter);
		}
	}
	
	delete reader;
	delete writer;

	return 0;
}
