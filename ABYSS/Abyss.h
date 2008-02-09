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

bool isCoordInternal(Coord4 c, Coord4 start, Coord4 size);

void outputSequences(const char* filename, PhaseSpace* pPS, Coord4 minCoord, Coord4 maxCoord);
void assemble(PhaseSpace* pPS, Coord4 minCoord, Coord4 maxCoord);
Sequence assembleSequence(PhaseSpace* pPS, 	PhaseSpaceBinIter sequenceIter);
Sequence BuildContig(PSequenceVector* extensions, PackedSeq& originalSeq);
void printUsage();


void trimSequences(PhaseSpace* pPS, Coord4 minCoord, Coord4 maxCoord);
void trimSequences2(PhaseSpace* pPS, Coord4 minCoord, Coord4 maxCoord, int trimNum);

#endif
