#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <fstream>
#include <string>
#include <math.h>
#include "PhaseSpace.h"
#include "CommonUtils.h"

typedef std::set<PackedSeq> PSeqSet;

int main(int argv, char** argc)
{
	if(argv < 3)
	{
		printf("Usage: estimateContigs <raw genome file> <kmer size>\n");
	}
	
	std::string rawSeqFile = argc[1];
	int kmerSize = atoi(argc[2]);

	
	// open file
	//printf("reading genome\n");
	std::ifstream handle;
	handle.open(rawSeqFile.c_str());

	handle.seekg (0, std::ios::end);
	int length = handle.tellg();
	handle.seekg (0, std::ios::beg);
	
	// allocate memory:
	char* buffer = new char[length];	
	
  	// read data as a block:
	handle.read (buffer,length);

	handle.close();	
	
	//printf("finished reading\n");
	
	Coord4 start;
	start.x = 0;
	start.y = 0;
	start.z = 0;
	start.w = 0;
	
	Coord4 size;
	size.x = kmerSize;
	size.y = kmerSize;
	size.z = kmerSize;
	size.w = kmerSize;
	
	PhaseSpace ps(kmerSize, start, size);
	
	//printf("starting phase space load\n");
	
	// make upper case
	for(int i = 0; i < length; i++)
	{
		buffer[i] = std::toupper(buffer[i]);
	}
	

	PSeqSet* seqSet = new PSeqSet;
	
	for(int i = 0; i < length - kmerSize; i++)
	{
		if(i % 100000 == 0)
		{
			//printf("loading %d/%d\n", i, length);
		}
		
		char* pStr = buffer + i;
		Sequence seq(pStr, kmerSize);
		//printf("seq: %s\n", seq.c_str());
		PackedSeq pSeq(seq);
		seqSet->insert(pSeq);
		//ps.addSequence(pSeq);
	}
	
	
	for(PSeqSet::iterator iter = seqSet->begin(); iter != seqSet->end(); iter++)
	{
		ps.addSequence(*iter);
	}
	
	delete seqSet;
	
	ps.finalizeBins();
	
	//printf("done load\n");
	
	int currSize = 0;
	for(int i = 0; i < length - kmerSize; i++)
	{
		char* pStr = buffer + i;
		Sequence seq(pStr, kmerSize);
		
		PackedSeq pSeq(seq);
		HitRecord hr = ps.calculateExtension(pSeq, SENSE);

		if(hr.getNumHits() == 1)
		{
			currSize++;
		}
		else
		{
			printf("%d\n", currSize);
			currSize = 0;
		}
	}
	delete [] buffer;
}
