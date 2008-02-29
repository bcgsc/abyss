#ifndef ASSEMBLYALGORITHMS_H
#define ASSEMBLYALGORITHMS_H

#include "ISequenceCollection.h"
#include "FastaReader.h"
#include "SeqRecord.h"

/*********************************************************
 * 
 * AssemblyAlgorithms.h
 * 
 * A collection of functions to operate on our sequence data
 * sets. These functions are designed to work with network
 * (parallel) or local (single cpu) data
 * 
 * 
 **********************************************************/
 
 
// Read a sequence file and load them into the collection
void loadSequences(ISequenceCollection* seqCollection, std::string fastaFile, int readLength, int kmerSize);

// Generate the adjacency information for all the sequences in the collection
// This is required before any other algorithm can run
void generateAdjacency(ISequenceCollection* seqCollection);

// trimming driver function, iteratively calls trimSequences to get rid of sequences that likely contain errors
void performTrim(ISequenceCollection* seqCollection, int readLen, int kmerSize);

// Function to perform the actual trimming. Walks the sequence space 
int trimSequences(ISequenceCollection* seqCollection, int maxBranchCull);

// The actual assembly function, takes in an ISequenceCollection pointer
void assemble(ISequenceCollection* seqCollection);

// BuildContig generates the sequence of a contig from the 2-element array of 
// sequence vectors (one for the antisense direction, one for sense direction)
Sequence BuildContig(PSequenceVector* extensions, const PackedSeq& originalSeq);

// remove the specified sequence from the collection and destroy its extensions
void removeSequenceAndExtensions(ISequenceCollection* seqCollection, const PackedSeq& seq);

#endif
