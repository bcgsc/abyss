#ifndef ABYSS_H
#define ABYSS_H

#include "Sequence.h"
#include "Reader.h"
#include "PathDriver.h"
#include "PairRecord.h"
#include "Writer.h"
#include "SeqRecord.h"
#include "Config.h"
#include "PartitionLoader.h"
#include "FastaWriter.h"
#include "FastaReader.h"
#include "PackedSeqWriter.h"
#include "SimpleSequenceSpace.h"


typedef std::pair<PackedSeq, extDirection> branchEnd;


bool isCoordInternal(Coord4 c, Coord4 start, Coord4 size);

void outputSequences(const char* filename, SimpleSequenceSpace* pSS, Coord4 minCoord, Coord4 maxCoord);
void assemble(SimpleSequenceSpace* pSS, Coord4 minCoord, Coord4 maxCoord);
void assemble2(SimpleSequenceSpace* pSS, Coord4 minCoord, Coord4 maxCoord);
Sequence assembleSequence(SimpleSequenceSpace* pSS, 	PhaseSpaceBinIter sequenceIter);
Sequence BuildContig(PSequenceVector& extensions, PackedSeq& originalSeq, extDirection dir);
Sequence BuildContig(PSequenceVector* extensions, PackedSeq& originalSeq);
void printUsage();


void trimSequences(SimpleSequenceSpace* pSS, Coord4 minCoord, Coord4 maxCoord);
void trimSequences2(SimpleSequenceSpace* pSS, Coord4 minCoord, Coord4 maxCoord, int maxBranchCull);
void outputBranchSizes(SimpleSequenceSpace* pSS, Coord4 minCoord, Coord4 maxCoord);

#endif
