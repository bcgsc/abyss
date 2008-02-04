#ifndef PATHWALKER_H
#define PATHWALKER_H

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

#endif