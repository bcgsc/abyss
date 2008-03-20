#ifndef COMMONDEFS_H
#define COMMONDEFS_H

#include <vector>
#include <set>
#include <map>
#include <stdlib.h>
#include <ext/hash_set>
#include <assert.h>
#include "ReadPrb.h"
#include "Prb.h"

const int MAX_FASTA_LINE = 1024;

enum FileMode
{
	FM_READ,
	FM_WRITE
	
};

enum FileType
{
	FT_FASTA, //ascii fasta file
	FT_SQB    //compressed sequence binary format
};

struct Coord4
{
	int x;
	int y;
	int z;
	int w;
};

// SENSE AND ANTISENSE HAVE TO BE ZERO AND ONE
enum extDirection
{
	SENSE = 0,
	ANTISENSE = 1,
	NUM_DIRECTIONS
};

enum CollectionState
{
	CS_LOADING,
	CS_FINALIZED
};

typedef std::vector<int> Count1D;
typedef std::vector<Count1D> Count2D;
typedef std::vector<Count2D> Count3D;
typedef std::vector<Count3D> Count4D;

const int NUM_BASES = 4;
const char BASES[NUM_BASES] = {'A', 'C', 'G', 'T'};

typedef std::string Sequence;

class PackedSeq;

// typedefs to make stl-based code somewhat more readable
typedef std::vector<Sequence>::const_iterator const_seq_iter;
typedef std::vector<Sequence>::iterator seq_iter;
typedef std::vector<Sequence>::const_reverse_iterator const_rev_seq_iter;
typedef std::vector<Sequence>::reverse_iterator rev_seq_iter;

typedef std::map<std::string, Sequence> SequenceMap;
typedef std::map<std::string, ReadPrb> PrbMap;


typedef std::vector<ReadPrb> PrbVector;

// The main data model for the program is a vector (contiguous array) of PackedSeqs
typedef std::vector<Sequence> SequenceVector;
typedef SequenceVector::const_iterator ConstSequenceVectorIterator;
typedef SequenceVector::iterator SequenceVectorIterator;

typedef std::vector<PackedSeq> PSequenceVector;
typedef PSequenceVector::const_iterator ConstPSequenceVectorIterator;
typedef PSequenceVector::iterator PSequenceVectorIterator;


struct PackedSeqEqual
{
	bool operator()(const PackedSeq& obj1, const PackedSeq& obj2) const;	
};

struct PackedSeqHasher
{
	size_t operator()(const PackedSeq& myObj) const;
};

typedef __gnu_cxx::hash_set<PackedSeq, PackedSeqHasher, PackedSeqEqual> SequenceDataHash;
//typedef std::set<PackedSeq> SequenceDataHash;
typedef SequenceDataHash::iterator SequenceCollectionHashIter;
typedef SequenceDataHash::const_iterator ConstSequenceCollectionHashIter;

typedef std::pair<SequenceCollectionHashIter, SequenceCollectionHashIter> SequenceHashIterPair;


typedef SequenceCollectionHashIter SequenceCollectionIterator;
//typedef PSequenceVectorIterator SequenceCollectionIterator;

#endif
