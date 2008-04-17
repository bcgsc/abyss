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

// Most operations are performed on the forward and reverse reads simulatenously, this structure holds the result of such operations
struct ResultPair
{
	bool forward;
	bool reverse;
};


struct Coord4
{
	int x;
	int y;
	int z;
	int w;
};

// IDs for reads and contigs
typedef int ReadID;
typedef std::string ContigID;

struct Position
{
	ContigID contig;
	int pos; // 0 indexed
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

// Definition of bases
const int NUM_BASES = 4;
const char BASES[NUM_BASES] = {'A', 'C', 'G', 'T'};

// A sequence is simply a c++ string
typedef std::string Sequence;
typedef std::vector<Sequence> SequenceVector;
typedef SequenceVector::const_iterator ConstSequenceVectorIterator;
typedef SequenceVector::iterator SequenceVectorIterator;
typedef std::map<std::string, Sequence> SequenceMap;

// Forward declare of a PackedSeq
class PackedSeq;
typedef std::vector<PackedSeq> PSequenceVector;
typedef std::set<PackedSeq> PSeqSet;
typedef PSequenceVector::const_iterator ConstPSequenceVectorIterator;
typedef PSequenceVector::iterator PSequenceVectorIterator;

// Contig definition
struct Contig
{
	Sequence seq;
	bool merged;
	bool repetitive;
};

// Contig typedefs
typedef std::map<ContigID, Contig> ContigMap;
typedef ContigMap::iterator CMIter;
typedef ContigMap::const_iterator ConstCMIter;


// Hash/Set functions
struct PackedSeqEqual
{
	bool operator()(const PackedSeq& obj1, const PackedSeq& obj2) const;	
};

struct PackedSeqHasher
{
	size_t operator()(const PackedSeq& myObj) const;
};


//
// Typedefs for the main data model
//
typedef __gnu_cxx::hash_set<PackedSeq, PackedSeqHasher, PackedSeqEqual> SequenceDataHash;
//typedef std::set<PackedSeq> SequenceDataHash;
typedef SequenceDataHash::iterator SequenceCollectionHashIter;
typedef SequenceDataHash::const_iterator ConstSequenceCollectionHashIter;

typedef std::pair<SequenceCollectionHashIter, SequenceCollectionHashIter> SequenceHashIterPair;


typedef SequenceCollectionHashIter SequenceCollectionIterator;
//typedef PSequenceVectorIterator SequenceCollectionIterator;

#endif
