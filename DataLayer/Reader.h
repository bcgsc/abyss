#ifndef READER_H
#define READER_H

#include <vector>
#include "CommonDefs.h"
#include "Sequence.h"
#include "PhaseSpace.h"
#include "ReadPrb.h"
#include "PackedSeq.h"

// Wrapper class to parse a variety of file types
class Reader
{
	
	public:
	
		Reader();
	
		// read in a fasta file and output the packed sequences into a vector 
		bool readFasta(const char* filename, PSequenceVector& outSequences) const;
		
		// read in an apb file and output the prb values into the map (id->prbvals)
		bool readAPB(const char* filename, PrbVector& outPrbs) const;
		
		// read in a phase space file
		bool readPhaseSpaceBinFile(const char* filename, PhaseSpace& phaseSpace) const;
	
};

#endif
