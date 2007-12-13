#ifndef COMMONDEFS_H
#define COMMONDEFS_H

#include <vector>
#include <map>
#include "ReadPrb.h"
#include "Prb.h"

enum FileMode
{
	FM_READ,
	FM_WRITE
	
};

typedef std::vector<int> Count1D;
typedef std::vector<Count1D> Count2D;
typedef std::vector<Count2D> Count3D;
typedef std::vector<Count3D> Count4D;

const int NUM_BASES = 4;
const char BASES[NUM_BASES] = {'A', 'C', 'G', 'T'};

typedef std::string Sequence;

typedef std::vector<Sequence>::const_iterator const_seq_iter;
typedef std::vector<Sequence>::iterator seq_iter;
typedef std::vector<Sequence>::const_reverse_iterator const_rev_seq_iter;
typedef std::vector<Sequence>::reverse_iterator rev_seq_iter;

typedef std::map<std::string, Sequence> SequenceMap;
typedef std::map<std::string, ReadPrb> PrbMap;

typedef std::vector<Sequence> SequenceVector;
typedef std::vector<ReadPrb> PrbVector;

typedef SequenceVector::const_iterator ConstSequenceVectorIterator;
typedef SequenceVector::iterator SequenceVectorIterator;

#endif
