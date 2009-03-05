#include "Aligner.h"
#include "PairedAlgorithms.h"
#include <cassert>
#include <cerrno>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <getopt.h>
#include <sstream>
#include <string>

using namespace std;

#define PROGRAM "KAligner"

static const char *VERSION_MESSAGE =
PROGRAM " (ABySS) " VERSION "\n"
"Written by Jared Simpson and Shaun Jackman.\n"
"\n"
"Copyright 2009 Canada's Michael Smith Genome Science Centre\n";

static const char *USAGE_MESSAGE =
"Usage: " PROGRAM " [OPTION]... QUERY TARGET\n"
"Align the sequences of QUERY against those of TARGET.\n"
"All perfect matches of at least k bases will be found.\n"
"\n"
"  -k, --kmer=KMER_SIZE  k-mer size\n"
"  -v, --verbose         display verbose output\n"
"      --help            display this help and exit\n"
"      --version         output version information and exit\n"
"\n"
"Report bugs to <" PACKAGE_BUGREPORT ">.\n";

namespace opt {
	static unsigned k;
	static int verbose;
}

static const char* shortopts = "k:o:v";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
	{ "kmer",        required_argument, NULL, 'k' },
	{ "verbose",     no_argument,       NULL, 'v' },
	{ "help",        no_argument,       NULL, OPT_HELP },
	{ "version",     no_argument,       NULL, OPT_VERSION },
	{ NULL, 0, NULL, 0 }
};


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
	bool die = false;
	for (char c; (c = getopt_long(argc, argv,
					shortopts, longopts, NULL)) != -1;) {
		istringstream arg(optarg != NULL ? optarg : "");
		switch (c) {
			case '?': die = true; break;
			case 'k': arg >> opt::k; break;
			case 'v': opt::verbose++; break;
			case OPT_HELP:
				cout << USAGE_MESSAGE;
				exit(EXIT_SUCCESS);
			case OPT_VERSION:
				cout << VERSION_MESSAGE;
				exit(EXIT_SUCCESS);
		}
	}

	if (opt::k <= 0) {
		cerr << PROGRAM ": missing -k,--kmer option\n";
		die = true;
	}

	if (argc - optind < 2) {
		cerr << PROGRAM ": missing arguments\n";
		die = true;
	} else if (argc - optind > 2) {
		cerr << PROGRAM ": too many arguments\n";
		die = true;
	}

	if (die) {
		cerr << "Try `" << PROGRAM
			<< " --help' for more information.\n";
		exit(EXIT_FAILURE);
	}

	string readsFile(argv[optind++]);
	string refFastaFile(argv[optind++]);

	if (opt::verbose > 0)
		cerr << "k: " << opt::k
			<< " Query: " << readsFile
			<< " Target: " << refFastaFile
			<< endl;

	Aligner aligner(opt::k);

	// Read the contigs into the ref DB
	readContigsIntoDB(refFastaFile, aligner);
	
	// Align the reads	
	alignReadsToDB(readsFile, aligner);
	
	return 0;
} 

static void assert_open(std::ifstream& f, const std::string& p)
{
	if (f.is_open())
		return;
	std::cerr << p << ": " << strerror(errno) << std::endl;
	exit(EXIT_FAILURE);
}

void readContigsIntoDB(std::string refFastaFile, Aligner& aligner)
{
	int count = 0;
	std::ifstream fileHandle(refFastaFile.c_str());	
	assert_open(fileHandle, refFastaFile);

	while(!fileHandle.eof() && fileHandle.peek() != EOF)
	{
		ContigID contigID;
		Sequence seq;
		int length;
		double coverage;
		
		PairedAlgorithms::parseContigFromFile(fileHandle, contigID, seq, length, coverage);
		aligner.addReferenceSequence(contigID, seq);

		count++;
		if (opt::verbose > 0
				&& count % 100000 == 0)
			cerr << "Read " << count << " contigs, "
				<< aligner.getNumSeqs() << " unique sequences\n";
	}
	if (opt::verbose > 0)
		cerr << "Read " << count << " contigs, "
			<< aligner.getNumSeqs() << " unique sequences\n";

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
	assert_open(fileHandle, readsFile);

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
