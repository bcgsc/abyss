#ifndef ASSEMBLYALGORITHMS_H
#define ASSEMBLYALGORITHMS_H

#include <set>
#include "ISequenceCollection.h"
#include "FastaReader.h"
#include "SeqRecord.h"
#include "FastaWriter.h"

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
struct Branch
{
	PSeqSet seqSet;
	PackedSeq lastSeq;
	
	void AddSequence(const PackedSeq& seq) { seqSet.insert(seq).first; lastSeq = seq;}
};


// Calculate the extensions for this sequence
HitRecord calculateExtension(ISequenceCollection* seqCollection, const PackedSeq& currSeq, extDirection dir);
 
// Read a sequence file and load them into the collection
void loadSequences(ISequenceCollection* seqCollection, std::string fastaFile, int readLength, int kmerSize);

// Generate the adjacency information for all the sequences in the collection
// This is required before any other algorithm can run
void generateAdjacency(ISequenceCollection* seqCollection);

// trimming driver function, iteratively calls trimSequences to get rid of sequences that likely contain errors
void performTrim(ISequenceCollection* seqCollection);

// Function to perform the actual trimming. Walks the sequence space 
int trimSequences(ISequenceCollection* seqCollection, int maxBranchCull);

// Pop bubbles (loops of sequence that diverge a single base, caused by SNPs or consistent sequence errors
void popBubbles(ISequenceCollection* seqCollection, int kmerSize);

// Remove extensions to/from ambiguous sequences to avoid generating redundant/wrong contigs
void splitAmbiguous(ISequenceCollection* seqCollection);

// The actual assembly function, takes in an ISequenceCollection pointer
void assemble(ISequenceCollection* seqCollection, int readLen, int kmerSize, IFileWriter* fileWriter);

// BuildContig generates the sequence of a contig from the 2-element array of 
// sequence vectors (one for the antisense direction, one for sense direction)
Sequence BuildContig(ISequenceCollection* seqCollection, HitVector* extensions, const PackedSeq& originalSeq, int contigID, int readLen, int kmerSize);

// remove the specified sequence from the collection and destroy its extensions
void removeSequenceAndExtensions(ISequenceCollection* seqCollection, const PackedSeq& seq);
void removeExtensionsToSequence(ISequenceCollection* seqCollection, const PackedSeq& seq, extDirection dir);

void outputSequences(const char* filename, ISequenceCollection* pSS);

#endif
