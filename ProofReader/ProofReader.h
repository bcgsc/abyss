#ifndef PROOFREADER_H
#define PROOFREADER_H

#include "Reader.h"
#include "SeqRecord.h"

// Correct the read set using the map of sequences to their multiplicity
// The corrected reads are output into the vector
void correctReads(const SeqRecord& seqMult, std::map<Sequence, Sequence>& corrections);

// Experimental correction that operates by filling holes where possible
void fillHoles(const SequenceVector& seqVector, const SeqRecord& multiplicity, const PhaseSpace& phase, std::map<Sequence, Sequence>& corrections);

// Correct a single sequence by examining the possibly permutations
Sequence correctSequence(const Sequence& seq, const SeqRecord& multiplicity);

// Output the reads
void outputCorrectedSequences(const SequenceVector& seqVector, const SeqRecord& multiplicity, std::map<Sequence, Sequence>& corrections);

#endif
