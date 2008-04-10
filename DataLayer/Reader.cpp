#include <fstream>
#include <iostream>
#include <stdio.h>

#include "Reader.h"

Reader::Reader()
{
	
}

bool Reader::readFasta(const char* filename, PSequenceVector& outSequences) const
{
	char headerBuffer[MAX_FASTA_LINE];
	char seqBuffer[MAX_FASTA_LINE];
	
	// open file and check that we can read from it
	std::ifstream inFile(filename);
	if(!inFile.is_open())
	{
		printf("Could not open fasta file %s\n", filename);
		return false;
	}
	
	char id[SEQUENCE_ID_LENGTH];
	
	while(!inFile.eof() && inFile.peek() != EOF)
	{
		// read in the header
		inFile.getline(headerBuffer, MAX_FASTA_LINE);
	
		// read in the sequence
		inFile.getline(seqBuffer, MAX_FASTA_LINE);
			
		// check if the reads were successful
		if(inFile.fail())
		{
			printf("error reading from fasta file %s. Is the line length too long?\n", filename);
			return false;
		}
		
		// parse the header
		if(sscanf(headerBuffer, ">%s %*s", id) != 1)
		{
			printf("invalid header format in %s, read failed\n", filename);
			return false;
		}
		
		// create the new sequence
		std::string idStr(id);
		Sequence newSeq(seqBuffer);
		PackedSeq pSeq(newSeq);
		outSequences.push_back(pSeq);
	}
	
	return true;
}

bool Reader::readPhaseSpaceBinFile(const char* filename, PhaseSpace& phaseSpace) const
{
	char headerBuffer[MAX_FASTA_LINE];
	char seqBuffer[MAX_FASTA_LINE];
	
	// open file and check that we can read from it
	std::ifstream inFile(filename);
	if(!inFile.is_open())
	{
		printf("Could not open fasta file %s\n", filename);
		return false;
	}
	
	while(!inFile.eof() && inFile.peek() != EOF)
	{
		// read in the header
		inFile.getline(headerBuffer, MAX_FASTA_LINE);
		
		// read in the sequence
		inFile.getline(seqBuffer, MAX_FASTA_LINE);
		
		// check if the reads were successful
		if(inFile.fail())
		{
			printf("error reading from fasta file %s. Is the line length too long?\n", filename);
			return false;
		}
		
		char id[SEQUENCE_ID_LENGTH];
		Coord4 c;
		// parse the header
		if(sscanf(headerBuffer, ">%s %d %d %d %d", id, &c.x, &c.y, &c.z, &c.w) != 5)
		{
			printf("invalid header format in %s, read failed\n", filename);
			return false;
		}

		// add the new sequence at the correct coordinate
		Sequence newSeq(seqBuffer);
		phaseSpace.addSequence(newSeq);
	}	
	
	return true;
}
