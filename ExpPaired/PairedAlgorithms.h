#ifndef PAIREDALGORITHMS_H
#define PAIREDALGORITHMS_H

#include "CommonDefs.h"
#include "PackedSeq.h"
#include "ISequenceCollection.h"
#include "AssemblyAlgorithms.h"
#include "ParentTree.h"
#include "PairRecord.h"
#include "PackedSeq.h"
#include "DirectedGraph.h"
#include "ContigData.h"

typedef DirectedGraph<ContigID, ContigData> ContigGraph;

namespace PairedAlgorithms
{

void readContigMap(std::string file, ContigMap& outMap);
void parseContigFromFile(std::ifstream& stream, ContigID& id, Sequence& seq, int& length, double& coverage);
void generateGraph(ContigGraph* pGraph, const ContigMap& contigMap, ISequenceCollection* pSC, size_t kmer, AlignmentCache* pDB);

};



#endif
