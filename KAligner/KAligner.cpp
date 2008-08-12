#include <stdio.h>
#include <math.h>
#include <iostream>
#include "SequenceCollectionHash.h"
#include "AssemblyAlgorithms.h"
#include "Options.h"
#include "FastaReader.h"
#include "Stats.h"
#include "PairUtils.h"
#include "PairedAlgorithms.h"
#include "Timer.h"

enum SequenceFormat
{
	SF_FASTA,
	SF_FASTQ
};

// Functions
void readContigsIntoDB(std::string refFastaFile, Aligner& aligner);
void alignReadsToDB(std::string readsFile, Aligner& aligner);

void readSequenceRecord(std::ifstream& stream, SequenceFormat type, std::string& readID, Sequence& seq);
void readFastaRecord(std::ifstream& stream, std::string& readID, Sequence& seq);
void readFastqRecord(std::ifstream& stream, std::string& readID, Sequence& seq);

int main(int argc, char** argv)
{
	if(argc < 4)
	{
		std::cout << "Usage: <kmer> <reads fasta | reads fastq> <ref fasta>\n";
		exit(1);
	}
	
	size_t kmer = atoi(argv[1]);
	std::string readsFile(argv[2]);
	std::string refFastaFile(argv[3]);

	std::cerr << "Kmer " << kmer << " Reads file: " << readsFile << " ref fasta: " << refFastaFile << std::endl;

	Aligner aligner(kmer);
	
	// Read the contigs into the ref DB
	readContigsIntoDB(refFastaFile, aligner);
	
	// Align the reads	
	alignReadsToDB(readsFile, aligner);
	
	return 0;
} 

void readContigsIntoDB(std::string refFastaFile, Aligner& aligner)
{
	int count = 0;
	std::ifstream fileHandle(refFastaFile.c_str());	
	while(!fileHandle.eof() && fileHandle.peek() != EOF)
	{
		ContigID contigID;
		Sequence seq;
		int length;
		double coverage;
		
		PairedAlgorithms::parseContigFromFile(fileHandle, contigID, seq, length, coverage);
		aligner.addReferenceSequence(contigID, seq);
		
		if(count % 100000 == 0)
		{
			std::cerr << "Read " << count << " contigs, " << aligner.getNumSeqs() << " seqs in the DB\n";
		}
		count++;
	}
	
	fileHandle.close();
}

void alignReadsToDB(std::string readsFile, Aligner& aligner)
{
	
	// Infer the reads file type
	SequenceFormat seqType;
	if(readsFile.find(".fq") != std::string::npos || readsFile.find(".fastq") != std::string::npos)
	{
		seqType = SF_FASTQ;
	}
	else if(readsFile.find(".fa") != std::string::npos || readsFile.find(".fasta") != std::string::npos)
	{
		seqType = SF_FASTA;
	}
	else
	{
		std::cerr << "Unknown file type!\n";
		assert(false);
	}
	
	std::ifstream fileHandle(readsFile.c_str());	
	while(!fileHandle.eof() && fileHandle.peek() != EOF)
	{
		std::string readID;
		Sequence readSeq;
		
		readSequenceRecord(fileHandle, seqType, readID, readSeq);
		
		// Filter for sequences only containing ACGT
		size_t pos = readSeq.find_first_not_of("ACGT");
		if (pos == std::string::npos) 
		{
			AlignmentVector avec;
			aligner.alignRead(readSeq, avec);
			if(!avec.empty())
			{
				std::cout << readID << "\t";
				for(AlignmentVector::iterator iter = avec.begin(); iter != avec.end(); ++iter)
				{
					std::cout << *iter << "\t";
				}
				std::cout << "\n";
			}		
		}
	}
}

void readSequenceRecord(std::ifstream& stream, SequenceFormat type, std::string& readID, Sequence& seq)
{
	switch(type)
	{
		case SF_FASTA:
			return readFastaRecord(stream, readID, seq);
		case SF_FASTQ:
			return readFastqRecord(stream, readID, seq);
		default:
			// unknown type;
			assert(false);
	}
}

void readFastaRecord(std::ifstream& stream, std::string& readID, Sequence& seq)
{
	std::string buffer;

	// Read the header.
	char c = stream.get();
	assert(c == '>');
	stream >> readID;
	
	// discard the remainder of the line
	getline(stream, buffer);
	
	// Read the sequence.
	getline(stream, seq);	
}

void readFastqRecord(std::ifstream& stream, std::string& readID, Sequence& seq)
{
	std::string buffer;

	// Read the header.
	char c = stream.get();
	assert(c == '@');
	stream >> readID;
	getline(stream, buffer);
	
	// Read the sequence.
	getline(stream, seq);

	// Read the quality values.
	getline(stream, buffer);
	assert(buffer[0] == '+');
	getline(stream, buffer);
}
