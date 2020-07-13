#ifndef RRESOLVER_BLOOMFILTERS_H
#define RRESOLVER_BLOOMFILTERS_H 1

#include "Common/Sequence.h"
#include "vendor/btl_bloomfilter/BloomFilter.hpp"

#include <vector>
#include <string>
#include <memory>
#include <cmath>
#include <map>

const int      HASH_NUM = 7; // For vanilla bloom filter
const int      SPACED_SEEDS_HASH_PER_SEED = 5;
const int      SPACED_SEEDS_COUNT = 6;
const int      SPACED_SEEDS_MIN_HITS = 1;
const double   SPACED_SEEDS_MISSES = 1;
const bool     SPACED_SEEDS_QC = true;
const double   SPACED_SEEDS_QC_MIN_BASES_PATTERN = 0.70;
const double   SPACED_SEEDS_QC_MIN_BASES_COMBINATION = 0.70;
const double   SPACED_SEEDS_QC_MIN_BASES_OVERALL = 0.65;
const double   SPACED_SEEDS_QC_MIN_ERRORS_OVERALL = 0.65;
const unsigned MAX_ERRORS_PER_RMER = 1;
const unsigned RMER_LOAD_STEP = 1;
const double   VANILLA_TO_SEEDS_MEM_RATIO = 1.15;

typedef std::vector<unsigned> SpacedSeed;
typedef std::vector<SpacedSeed> SpacedSeeds;

class VanillaBloomFilter : protected BloomFilter {

public:

	VanillaBloomFilter(size_t filterSize, int kmerSize): BloomFilter(filterSize, HASH_NUM, kmerSize) {}

  size_t getPop() const;
  double getFPR() const { return std::pow(double(getPop())/double(m_size), double(m_hashNum)); }

  unsigned getHashNum() const { return BloomFilter::getHashNum(); }
	unsigned getKmerSize() const { return BloomFilter::getKmerSize(); }
  size_t getFilterSize() const { return BloomFilter::getFilterSize(); }

  void loadSequence(const Sequence& sequence);

  bool contains(const size_t precomputed[]) const { return BloomFilter::contains(precomputed); }
  bool contains(vector<size_t> const &precomputed) const { return BloomFilter::contains(precomputed); }

};

class SpacedSeedsBloomFilter : protected BloomFilter {

public:

	SpacedSeedsBloomFilter(size_t filterSize, int kmerSize);

  size_t getPop() const;
  double getFPR() const;

  unsigned getHashNum() const { return BloomFilter::getHashNum(); }
	unsigned getKmerSize() const { return BloomFilter::getKmerSize(); }
  size_t getFilterSize() const { return BloomFilter::getFilterSize(); }

  void loadSequence(const Sequence& sequence);

  SpacedSeeds getHitSeeds(const size_t precomputed[]) const {
    SpacedSeeds hitSeeds;
    size_t normHash;
    unsigned seedHits;
    for (unsigned i = 0; i < SPACED_SEEDS_COUNT; i++) {
      seedHits = 0;
      for (unsigned j = 0; j < SPACED_SEEDS_HASH_PER_SEED; j++) {
        normHash = precomputed[i * SPACED_SEEDS_HASH_PER_SEED + j] % m_size;
        if (m_filter[normHash / bitsPerChar] & bitMask[normHash % bitsPerChar]) { seedHits++; }
      }
      if (seedHits == SPACED_SEEDS_HASH_PER_SEED) {
        hitSeeds.push_back(spacedSeeds[i]);
      }
    }
    return hitSeeds;
  }

  SpacedSeeds getSpacedSeeds() const { return spacedSeeds; }

private:

  SpacedSeeds spacedSeeds;

};

extern VanillaBloomFilter* g_vanillaBloom;
extern SpacedSeedsBloomFilter* g_spacedSeedsBloom;

void buildFilters(const std::vector<std::string>& readFilepaths, const int r, const size_t bloomBytesTotal);

#endif
