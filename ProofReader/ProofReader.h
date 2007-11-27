#ifndef PROOFREADER_H
#define PROOFREADER_H

#include "Reader.h"
#include "SeqRecord.h"

void correctReads(SequenceVector& seqVector, const PrbVector& prbVector, const SeqRecord& multiplicity);
void fillHoles(const SequenceVector& seqVector, const PrbVector& prbVector, const SeqRecord& multiplicity, std::map<Sequence, Sequence>& corrections);

Sequence correctSequence(const Sequence& seq, const ReadPrb& readPrb, const SeqRecord& multiplicity);
void outputReadSet(const SequenceVector& seqVector, const PrbVector& prbVector, const SeqRecord& multiplicity, std::map<Sequence, Sequence>& corrections);

#endif
