#ifndef ASSEMBLY_SEQUENCECOLLECTION_H
#define ASSEMBLY_SEQUENCECOLLECTION_H 1

#include "config.h"
#include "Kmer.h"
#include "KmerData.h"

#if HAVE_GOOGLE_SPARSE_HASH_MAP
# include <google/sparse_hash_map>
typedef google::sparse_hash_map<Kmer, KmerData, hash<Kmer> >
	SequenceDataHash;
#else
# include "UnorderedMap.h"
typedef unordered_map<Kmer, KmerData, hash<Kmer> >
	SequenceDataHash;
#endif

#include "Assembly/ISequenceCollection.h"
#include "Assembly/DBG.h"

#endif
