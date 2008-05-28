#ifndef PAIREDALGORITHMS_H
#define PAIREDALGORITHMS_H

#include "CommonDefs.h"
#include "PackedSeq.h"
#include "ISequenceCollection.h"
#include "AssemblyAlgorithms.h"
#include "ParentTree.h"
#include "PairRecord.h"
#include "PackedSeq.h"

enum ContigStartState
{
	CSS_ISLAND, // open on both ends
	CSS_ENDPOINT, // one end open
	CSS_AMBIGUOUS // ambiguous	
};

struct ContigStart
{
	PackedSeq seq;
	ContigStartState state;
	extDirection dir; // only valid for state = SC_ENDPOINT
};

typedef std::vector<ContigStart> contigStartVec;

namespace PairedAlgorithms
{

void ReadContigs(std::string file, ContigMap& outMap);
void generateStartList(ISequenceCollection* seqCollection, contigStartVec& startList);
bool processNonlinearExtensionForBranch(ISequenceCollection* seqCollection, PairRecord* pPairRecord, BranchRecord& branch, PackedSeq& currSeq, ExtensionRecord extensions, int multiplicity);
bool deconvolvePaths(ISequenceCollection* seqCollection, PairRecord* pPairRecord, const BranchRecord& currBranch, const PackedSeq& branchPoint, extDirection dir, PackedSeq& chosenNode);

};



#endif
