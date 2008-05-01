#ifndef ASSEMBLYALGORITHMS_H
#define ASSEMBLYALGORITHMS_H

#include <set>
#include "ISequenceCollection.h"
#include "FastaReader.h"
#include "SeqRecord.h"
#include "FastaWriter.h"
#include "BranchGroup.h"
#include "BranchRecord.h"

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
enum TrimStatus
{
	TS_NOTRIM, // sequence is contiguious on both ends
	TS_ISLAND, // sequence is completely isolated
	TS_TRIMMABLE // sequence can be trimmed
};

//
//
// Data preperation functions
//
//

// Read a sequence file and load them into the collection
void loadSequences(ISequenceCollection* seqCollection, std::string fastaFile, int readLength, int kmerSize);

// Generate the adjacency information for all the sequences in the collection
// This is required before any other algorithm can run
void generateAdjacency(ISequenceCollection* seqCollection);

//
//
// Trimming (error removal functions)
//
//

// trimming driver function, iteratively calls trimSequences to get rid of sequences that likely contain errors
void performTrim(ISequenceCollection* seqCollection);

// Function to perform the actual trimming. Walks the sequence space 
int trimSequences(ISequenceCollection* seqCollection, int maxBranchCull);

// Check whether a sequence can be trimmed
TrimStatus checkSeqForTrim(ISequenceCollection* seqCollection, const PackedSeq& seq, extDirection& outDir);

// process a terminated branch for trimming
bool processTerminatedBranch(ISequenceCollection* seqCollection, BranchRecord& branch);

// Process the extensions of the current sequence for trimming
bool processExtensionForBranchTrim(BranchRecord& branch, PackedSeq& currSeq, ExtensionRecord extensions);

// Polymorphism removal

// Pop bubbles (loops of sequence that diverge a single base, caused by SNPs or consistent sequence errors
int popBubbles(ISequenceCollection* seqCollection, int kmerSize);

// Populate the branch group with the initial extensions to this sequence
void initiateBranchGroup(BranchGroup& group, const PackedSeq& seq, const SeqExt& extension, size_t maxBubbleSize);

// process an a branch group extension
void processBranchGroupExtension(BranchGroup& group, size_t branchIndex, const PackedSeq& seq, ExtensionRecord extensions);

// collapse bubbles that are joined together
void collapseJoinedBranches(ISequenceCollection* seqCollection, BranchGroup& group);

//
//
// Split the remaining ambiguous nodes to allow for a non-redundant assembly
//
//

// Remove extensions to/from ambiguous sequences to avoid generating redundant/wrong contigs
void splitAmbiguous(ISequenceCollection* seqCollection);

//
//
// Assembly functions
// 
//

// The actual assembly function, takes in an ISequenceCollection pointer
void assemble(ISequenceCollection* seqCollection, int readLen, int kmerSize, IFileWriter* fileWriter);

// BuildContig generates the sequence of a contig from the 2-element array of 
// sequence vectors (one for the antisense direction, one for sense direction)
Sequence BuildContig(ISequenceCollection* seqCollection, HitVector* extensions, const PackedSeq& originalSeq, int contigID, int readLen, int kmerSize);

//
//
// Generic functions to operate on the data set
//
//
// remove the specified sequence from the collection and destroy its extensions
void removeSequenceAndExtensions(ISequenceCollection* seqCollection, const PackedSeq& seq);
void removeExtensionsToSequence(ISequenceCollection* seqCollection, const PackedSeq& seq, extDirection dir);

// Calculate the extensions for this sequence
HitRecord calculateExtension(ISequenceCollection* seqCollection, const PackedSeq& currSeq, extDirection dir);

// Generate all the sequences for the extension record
void generateSequencesFromExtension(const PackedSeq& currSeq, extDirection dir, SeqExt extension, PSequenceVector& outseqs);

// Output all the sequences that remain (have not been deleted) in the dataset
void outputSequences(const char* filename, ISequenceCollection* pSS);

#endif
