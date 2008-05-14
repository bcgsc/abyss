#ifndef ABYSS_H
#define ABYSS_H

#include "Sequence.h"
#include "Reader.h"
#include "PairRecord.h"
#include "Writer.h"
#include "SeqRecord.h"
#include "Config.h"
#include "FastaWriter.h"
#include "FastaReader.h"
#include "PackedSeqWriter.h"
#include "AssemblyData.h"
#include "BranchRecord.h"

struct PairScore
{
	double weight;
	PackedSeq seq;	
};

struct ContigStart
{
	PackedSeq seq;
	extDirection dir;
};

typedef std::list<ContigStart> startList;

typedef std::vector<PairScore> PairScoreVec;
				
void loadPairs(PairRecord& pairRec, std::string pairFile, int kmerSize);
int assemble_paired(ISequenceCollection* seqCollection, PairRecord& pairRecord, FastaWriter* originalWriter, FastaWriter* newWriter);
int assemble2(ISequenceCollection* seqCollection, PairRecord& pairRecord, FastaWriter* originalWriter, FastaWriter* newWriter);

BranchState processExtensionPaired(BranchRecord& branch, PackedSeq& currSeq, ExtensionRecord extensions, bool inBranchDetected);
void loadSeqHash(std::string packedSeqFile, ISequenceCollection* pSS);

#endif
