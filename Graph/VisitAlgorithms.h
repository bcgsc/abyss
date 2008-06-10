#ifndef VISITALGORITHMS_H
#define VISITALGORITHMS_H

#include "CommonDefs.h"
#include "FastaWriter.h"
#include "PackedSeq.h"
#include "PairRecord.h"
#include "ContigData.h"
#include "AlignmentCache.h"
#include "DirectedGraph.h"

struct ContigDataFunctions
{
	ContigDataFunctions(size_t k ,AlignmentCache* pAliDB) : m_kmer(k), m_pDatabase(pAliDB) { m_overlap = m_kmer - 1; }
	
	
	void merge(const ContigID& targetKey, ContigData& targetData, const ContigID& slaveKey, const ContigData& slaveData, extDirection dir, bool reverse, bool removeMerge);
	void deleteCallback(const ContigID& slaveKey, const ContigData& slaveData);
	
	bool check(ContigData& target, const ContigData& slave, extDirection dir, bool reverse);
	
	size_t m_kmer;
	AlignmentCache* m_pDatabase;
	size_t m_overlap;
	
};

struct DBGenerator
{
	DBGenerator(AlignmentCache* pDB) : m_pDatabase(pDB) { }
	void visit(const ContigID& id, const ContigData& data);
	
	AlignmentCache* m_pDatabase;
	
	
};


struct ContigDataOutputter
{
	ContigDataOutputter(FastaWriter* pWriter) : m_pWriter(pWriter), m_numOutput(0) { }
	void visit(const ContigID& id, const ContigData& data);
	FastaWriter* m_pWriter;
	size_t m_numOutput;
};

struct PairedResolver
{
	PairedResolver(size_t kmer, size_t maxLength, PairRecord* pPairs, AlignmentCache* pDB) : m_kmer(kmer), m_maxlength(maxLength), m_pPairs(pPairs), m_pDatabase(pDB) { }
		
	int resolve(DirectedGraph<ContigID, ContigData>* pGraph, const ContigID id);
	
	size_t scorePairsToSet(PSequenceVector& pairs, PSeqSet& seqSet);
	
	size_t m_kmer;
	size_t m_maxlength;
	PairRecord* m_pPairs;
	AlignmentCache* m_pDatabase;
	
};

struct LinkFinder
{
	LinkFinder(PairRecord* pPairs, AlignmentCache* pDB) : m_pPairs(pPairs), m_pDatabase(pDB) { }
	
	void visit(const ContigID& id, ContigData& data);
	void getReachableSet(const ContigID& id, const ContigData& data, extDirection dir, ContigIDSet& idSet);
	
	PairRecord* m_pPairs;
	AlignmentCache* m_pDatabase;
	
};

struct SequenceDataCost
{
	SequenceDataCost(size_t kmer) : m_overlap(kmer - 1) { }
	
	size_t cost(const ContigData& data) { return data.m_seq.length() - m_overlap; } 
	size_t m_overlap;
};

#endif
