#ifndef ASSEMBLYALGORITHMS_H
#define ASSEMBLYALGORITHMS_H

#include "BranchGroup.h"
#include "BranchRecord.h"
#include "FastaWriter.h"
#include "ISequenceCollection.h"

/*********************************************************
 * 
 * AssemblyAlgorithms.h
 * 
 * A collection of functions to operate on sequence data
 * sets. These functions are designed to work with network
 * (parallel) or local (single cpu) data
 * 
 * 
 **********************************************************/
enum SeqContiguity
{
	SC_INVALID, // sequence has been deleted/seen
	SC_ISLAND, // sequence is completely isolated
	SC_ENDPOINT, // one end of the sequence is open 
	SC_CONTIGUOUS // the sequence is closed on both ends
};



namespace AssemblyAlgorithms
{

//
//
// Data preperation functions
//
//

// Read a sequence file and load them into the collection
void loadSequences(ISequenceCollection* seqCollection,
		std::string inFile);

// Generate the adjacency information for all the sequences in the collection
// This is required before any other algorithm can run
void generateAdjacency(ISequenceCollection* seqCollection);

/** Return the first local minimum of the k-mer coverage histogram. */
unsigned minimumCoverage(/*const*/ ISequenceCollection& c);

/* Erosion. Remove k-mer from the ends of blunt contigs. */
unsigned erodeEnds(ISequenceCollection* seqCollection);
unsigned erode(ISequenceCollection* c, const PackedSeq& seq);
unsigned getNumEroded();

// trimming driver function, iteratively calls trimSequences to get rid of sequences that likely contain errors
void performTrim(ISequenceCollection* seqCollection, int start = 1);

// Function to perform the actual trimming. Walks the sequence space 
int trimSequences(ISequenceCollection* seqCollection, int maxBranchCull);
unsigned removeMarked(ISequenceCollection* pSC);

// Check whether a sequence can be trimmed
SeqContiguity checkSeqContiguity(ISequenceCollection* seqCollection, const PackedSeq& seq, extDirection& outDir);

// process a terminated branch for trimming
bool processTerminatedBranchTrim(ISequenceCollection* seqCollection, BranchRecord& branch);

// Process the extensions of the current sequence for trimming
bool processLinearExtensionForBranch(BranchRecord& branch, PackedSeq& currSeq, ExtensionRecord extensions, int multiplicity);

// Polymorphism removal

// Pop bubbles (loops of sequence that diverge a single base, caused by SNPs or consistent sequence errors
int popBubbles(ISequenceCollection* seqCollection, int kmerSize);

// Populate the branch group with the initial extensions to this sequence
void initiateBranchGroup(BranchGroup& group, const PackedSeq& seq, const SeqExt& extension, int multiplicity, size_t maxBubbleSize);

// process an a branch group extension
bool processBranchGroupExtension(BranchGroup& group, size_t branchIndex, const PackedSeq& seq, ExtensionRecord extensions, int multiplicity);

/** Write SNP information to a file. */
void writeSNP(BranchGroup& group, unsigned id);

// collapse bubbles that are joined together
void collapseJoinedBranches(ISequenceCollection* seqCollection, BranchGroup& group);

//
//
// Split the remaining ambiguous nodes to allow for a non-redundant assembly
//
//

// Remove extensions to/from ambiguous sequences to avoid generating redundant/wrong contigs
unsigned markAmbiguous(ISequenceCollection* seqCollection);
unsigned splitAmbiguous(ISequenceCollection* seqCollection);

unsigned assembleContig(ISequenceCollection* seqCollection,
		FastaWriter* writer, BranchRecord& branch, unsigned id);
unsigned assemble(ISequenceCollection* seqCollection,
		FastaWriter* fileWriter = NULL);

// remove the specified sequence from the collection and destroy its extensions
void removeSequenceAndExtensions(ISequenceCollection* seqCollection, const PackedSeq& seq);
void removeExtensionsToSequence(ISequenceCollection* seqCollection, const PackedSeq& seq, extDirection dir);

// Generate all the sequences for the extension record
void generateSequencesFromExtension(const PackedSeq& currSeq, extDirection dir, SeqExt extension, PSequenceVector& outseqs);

};

#endif
