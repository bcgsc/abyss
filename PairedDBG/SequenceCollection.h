#ifndef PAIREDDBG_SEQUENCECOLLECTION_H
#define PAIREDDBG_SEQUENCECOLLECTION_H 1

#include "config.h"
#include "Dinuc.h"
#include "KmerPair.h"
#include "Assembly/VertexData.h"

typedef VertexData<Dinuc, DinucSet> KmerPairData;

#if HAVE_GOOGLE_SPARSE_HASH_MAP
# include <google/sparse_hash_map>
typedef google::sparse_hash_map<KmerPair, KmerPairData, hash<KmerPair> >
	SequenceDataHash;
#else
# include "Common/UnorderedMap.h"
typedef unordered_map<KmerPair, KmerPairData, hash<KmerPair> >
	SequenceDataHash;
#endif

#include "Assembly/DBG.h"
#include "PairedDBG/BranchRecord.h"

#endif
