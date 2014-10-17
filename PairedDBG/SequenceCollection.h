#ifndef PAIREDDBG_SEQUENCECOLLECTION_H
#define PAIREDDBG_SEQUENCECOLLECTION_H 1

#include "config.h"
#include "KmerPair.h"
#include "KmerPairData.h"

#if HAVE_GOOGLE_SPARSE_HASH_MAP
# include <google/sparse_hash_map>
typedef google::sparse_hash_map<KmerPair, KmerPairData, hash<KmerPair> >
	SequenceDataHash;
#else
# include "UnorderedMap.h"
typedef unordered_map<KmerPair, KmerPairData, hash<KmerPair> >
	SequenceDataHash;
#endif

#include "Assembly/ISequenceCollection.h"
#include "Assembly/DBG.h"

#endif
