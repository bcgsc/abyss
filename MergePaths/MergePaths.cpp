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

typedef std::vector<MergeNode> MergeNodeVector;

struct PathMergeRecord
{
	ContigID rootID;
	extDirection dir;
	
	MergeNodeVector mergeNodeVec;
	
};

// Functions
void readIDIntPair(std::string str, ContigID& id, int& i);
void mergePaths(std::string pathFile, ContigMap sourceContigs, size_t kmer);
void parsePathLine(std::string pathLine, PathMergeRecord& mergeRecord);
void mergeSequences(Sequence& rootContig, const Sequence& otherContig, extDirection dir, bool isReversed, size_t kmer);

int main(int argc, char** argv)
{
	if(argc < 4)
	{
		std::cout << "Usage: <kmer> <contigFile> <paths file>\n";
		exit(1);
	}
	
	size_t kmer = atoi(argv[1]);
	std::string contigFile(argv[2]);
	std::string pathFile(argv[3]);

	std::cout << "Kmer " << kmer << " Contig File: " << contigFile << " path file: " << pathFile << std::endl;

	// Read the contigs
	
	// Set up the ID->sequence mapping
	ContigMap contigMap;
	PairedAlgorithms::readContigMap(contigFile, contigMap);
	
	// Read the paths file and merge as necessary
	mergePaths(pathFile, contigMap, kmer);
	return 1;
} 

void mergePaths(std::string pathFile, ContigMap sourceContigs, size_t kmer)
{
	std::ifstream pathStream(pathFile.c_str());
	FastaWriter writer("final.fa");
	int count = 0;
	
	while(!pathStream.eof() && pathStream.peek() != EOF)
	{
		// read a line
		std::string pathRecord;
		getline(pathStream, pathRecord);
		
		// parse the line
		PathMergeRecord mergeRec;
		parsePathLine(pathRecord, mergeRec);
		
		std::cout << "Attempting to merge " << mergeRec.rootID << " in dir " << mergeRec.dir << "\n";
		Sequence merged = sourceContigs[mergeRec.rootID].seq;
		assert(!merged.empty());
		
		bool flip = false;
		for(MergeNodeVector::iterator iter = mergeRec.mergeNodeVec.begin(); iter != mergeRec.mergeNodeVec.end(); ++iter)
		{
			flip = flip ^ iter->isRC;
			std::cout << "	merging in " << iter->id << "(" << flip << ")\n";
			Sequence otherContig = sourceContigs[iter->id].seq;
			assert(!otherContig.empty());
			mergeSequences(merged, otherContig, mergeRec.dir, flip, kmer);
		}
		
		writer.WriteSequence(merged, count, 0.0f);
		count++;
		
	}
	
	pathStream.close();
}


//
//
//
void mergeSequences(Sequence& rootContig, const Sequence& otherContig, extDirection dir, bool isReversed, size_t kmer)
{
	size_t overlap = kmer - 1;
	
	// should the slave be reversed?
	Sequence slaveSeq = otherContig;
	if(isReversed)
	{
		slaveSeq = reverseComplement(slaveSeq);
	}
	
	const Sequence* leftSeq;
	const Sequence* rightSeq;
	// Order the contigs
	if(dir == SENSE)
	{
		leftSeq = &rootContig;
		rightSeq = &slaveSeq;
	}
	else
	{
		leftSeq = &slaveSeq;
		rightSeq = &rootContig;
	}
	
	// Get the last k bases of the left and the first k bases of the right
	PackedSeq leftEnd = leftSeq->substr(leftSeq->length() - overlap, overlap);
	PackedSeq rightBegin = rightSeq->substr(0, overlap);

	// ensure that there is a legitimate k-1 overlap between these sequences	
	if(leftEnd != rightBegin)
	{
		printf("merge called data1: %s %s (%d, %d)\n", rootContig.c_str(), otherContig.c_str(), dir, isReversed);	
		printf("left end %s, right begin %s\n", leftEnd.decode().c_str(), rightBegin.decode().c_str());
		assert(leftEnd == rightBegin);
	}
	
	// TODO: make this in-place?
	// generate the merged sequence
	Sequence merged = *leftSeq;
	merged.append(rightSeq->substr(overlap));
	rootContig = merged;
}

void parsePathLine(std::string pathLine, PathMergeRecord& mergeRecord)
{
	std::string discard;
	
	// discard the seperator
	std::stringstream pStream(pathLine);	
	pStream >> discard;
	
	// read in the root info
	std::string rootInfo;
	pStream >> rootInfo;
	
	ContigID rootID;
	int rootDir;
	readIDIntPair(rootInfo, rootID, rootDir);
	
	mergeRecord.rootID = rootID;
	mergeRecord.dir = (extDirection)rootDir;
	
	// discard the seperator
	pStream >> discard;
	
	// read all the fields
	bool stop = false;
	while(!stop)
	{
		// Read a record
		std::string record;
		pStream >> record;
		
		if(record.empty())
		{
			stop = true;
		}
		else
		{
			// parse the record
			ContigID id;
			int isRC;
			readIDIntPair(record, id, isRC);
			
			MergeNode node;
			node.id = id;
			node.isRC = (bool)isRC;
			
			mergeRecord.mergeNodeVec.push_back(node);
		}
	}
}


void readIDIntPair(std::string str, ContigID& id, int& i)
{
	std::stringstream ss(str);
	
	// read in the id
	getline(ss, id, ',');
	
	std::string intStr;
	ss >> intStr;
	
	std::stringstream convertor(intStr);
	convertor >> i;
}

