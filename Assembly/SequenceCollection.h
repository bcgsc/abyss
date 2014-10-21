#ifndef ASSEMBLY_SEQUENCECOLLECTION_H
#define ASSEMBLY_SEQUENCECOLLECTION_H 1

#include "config.h"
#include "SeqExt.h"
#include "VertexData.h"
#include "Common/Kmer.h"

typedef VertexData<uint8_t, SeqExt> KmerData;
typedef KmerData::SymbolSetPair ExtensionRecord;

#if HAVE_GOOGLE_SPARSE_HASH_MAP
# include <google/sparse_hash_map>
typedef google::sparse_hash_map<Kmer, KmerData, hash<Kmer> >
	SequenceDataHash;
#else
# include "Common/UnorderedMap.h"
typedef unordered_map<Kmer, KmerData, hash<Kmer> >
	SequenceDataHash;
#endif

#include "Assembly/DBG.h"
#include "Assembly/BranchRecord.h"

#endif
