#ifndef RRESOLVER_BLOOMFILTERS_H
#define RRESOLVER_BLOOMFILTERS_H 1

#include "btllib/bloom_filter.hpp"

#include <cmath>
#include <map>
#include <memory>
#include <string>
#include <vector>

const int HASH_NUM = 7; // For vanilla bloom filter
const int SPACED_SEEDS_HASHES_PER_SEED = 5;
const int SPACED_SEEDS_COUNT = 6;
const int SPACED_SEEDS_MIN_HITS = 1;
const double SPACED_SEEDS_MISSES = 1;
const bool SPACED_SEEDS_QC = true;
const double SPACED_SEEDS_QC_MIN_BASES_PATTERN = 0.70;
const double SPACED_SEEDS_QC_MIN_BASES_COMBINATION = 0.70;
const double SPACED_SEEDS_QC_MIN_BASES_OVERALL = 0.65;
const double SPACED_SEEDS_QC_MIN_ERRORS_OVERALL = 0.65;
const unsigned MAX_ERRORS_PER_RMER = 1;
const double VANILLA_TO_SEEDS_MEM_RATIO = 1.15;

extern btllib::KmerBloomFilter *g_vanillaBloom;
extern btllib::SeedBloomFilter *g_spacedSeedsBloom;

void
buildFilters(
    const std::vector<std::string>& read_filepaths,
    const int r,
    const size_t bloom_bytes_total);

#endif
