#include "BloomFilters.h"

#include "RUtils.h"
#include "RAlgorithmsShort.h"
#include "Common/IOUtil.h"
#include "Common/Options.h"
#include "DataLayer/FastaReader.h"
#include "vendor/nthash/stHashIterator.hpp"
#include "vendor/nthash/ntHashIterator.hpp"

#include "btllib/include/btllib/seq_reader.hpp"

#include <algorithm>
#include <random>
#include <cassert>
#include <cmath>

static std::vector<std::string>
generateSpacedSeedsPatterns(const int count, const int size, const int misses) {
  assert(count < size);
  assert(misses < count);
  assert(misses > 0);

  std::vector<std::string> seeds;
  auto randomEngine = std::mt19937();
  std::vector<int> seedsPermutation;

  for (int i = 0; i < count; i++) {
    seeds.push_back(std::string(size, '1'));
    seedsPermutation.push_back(i);
  }

  for (int i = 0; i < (size + 1) / 2; i++) {
    std::shuffle(seedsPermutation.begin(), seedsPermutation.end(), randomEngine);
    for (int j = 0; j < count; j++) {
      char c = j < misses ? '0' : '1';
      seeds[seedsPermutation[j]][i] = c;
    }
    if (i < size / 2) {
      for (int j = 0; j < count; j++) {
        seeds[count - j - 1][size - i - 1] = seeds[j][i];
      }
    }
  }

  return seeds;
}

void QCSpacedSeedsPatterns(const std::vector<std::string> &patterns) {
  int k = patterns[0].size();

  for (const auto &pattern : patterns) {
    int patternBasesCovered = k;
    bool hasZero = false;
    for (char c : pattern) {
      if (c == '0') {
        hasZero = true;
        patternBasesCovered--;
        break;
      }
    }
    if (opt::verbose && !hasZero) {
      std::cerr << "A spaced seed has no zeros!\n";
    }
    if (opt::verbose && (patternBasesCovered < k * SPACED_SEEDS_QC_MIN_BASES_PATTERN)) {
      std::cerr << "A spaced seed pattern does not cover enough bases!\n";
    }
  }

  auto combinations = genCombinations(SPACED_SEEDS_COUNT, SPACED_SEEDS_MIN_HITS);
  std::string overallErrorCoverage(k, '1');
  int overallErrorsCovered = 0;
  std::string overallBaseCoverage(k, '0');
  int overallBasesCovered = 0;
  std::string worstCombinationCoverage(k, '1');
  int worstCombinationBasesCovered = k;

  for (auto combination : combinations) {
    std::string combinationCoverage(k, '0');
    int combinationBasesCovered = 0;
    for (auto index : combination) {
      const auto &pattern = patterns[index];
      for (unsigned i = 0; i < pattern.size(); i++) {
        if (pattern[i] == '1' && combinationCoverage[i] != '1') {
          combinationCoverage[i] = '1';
          combinationBasesCovered++;
        }
      }
    }
    if (combinationBasesCovered < worstCombinationBasesCovered) {
      worstCombinationCoverage = combinationCoverage;
      worstCombinationBasesCovered = combinationBasesCovered;
    }
    if (opt::verbose && (combinationBasesCovered < combinationCoverage.size() * SPACED_SEEDS_QC_MIN_BASES_COMBINATION)) {
      std::cerr << "A spaced seed combination does not cover enough bases!\n";
    }

    for (unsigned i = 0; i < combinationCoverage.size(); i++) {
      if (combinationCoverage[i] == '0') {
        if (overallErrorCoverage[i] != '0') {
          overallErrorCoverage[i] = '0';
          overallErrorsCovered++;
        }
      } else {
        if (overallBaseCoverage[i] != '1') {
          overallBaseCoverage[i] = '1';
          overallBasesCovered++;
        }
      }
    }
  }

  if (opt::verbose) {
    std::cerr << std::fixed;
    std::cerr << "Worst combination coverage: " << worstCombinationCoverage << '\n';
    std::cerr << "Worst combination bases covered: " << worstCombinationBasesCovered / double(k) * 100.0 << "%\n";
    std::cerr << "Overall base coverage:\n" << overallBaseCoverage << '\n';
    std::cerr << "Bases covered: " << overallBasesCovered / double(k) * 100.0 << "%\n";
    std::cerr << "Overall error coverage:\n" << overallErrorCoverage << '\n';
    std::cerr << "Errors covered: " << overallErrorsCovered / double(k) * 100.0 << "%\n";
    std::cerr << std::defaultfloat << std::endl;
  }

  if (opt::verbose && (overallErrorsCovered < k * SPACED_SEEDS_QC_MIN_ERRORS_OVERALL)) {
    std::cerr << "Spaced seeds do not cover enough error positions!\n\n";
  }
  if (opt::verbose && (overallBasesCovered < k * SPACED_SEEDS_QC_MIN_BASES_OVERALL)) {
    std::cerr << "Spaced seeds do not cover enough base positions!\n\n";
  }
  std::cerr << std::flush;
}

size_t VanillaBloomFilter::getPop() const {
  size_t i, popBF = 0;
  #pragma omp parallel for reduction(+:popBF)
  for (i = 0; i < (m_size + 7) / 8; i++) {
    popBF = popBF + popCnt(m_filter[i]);
  }
  return popBF;
}

void VanillaBloomFilter::loadSequence(const Sequence& sequence) {
  if (sequence.size() >= m_kmerSize) {
#if RMER_LOAD_STEP > 1
    unsigned offset = 0;
#endif
    for (ntHashIterator it(sequence, HASH_NUM, m_kmerSize); it != ntHashIterator::end();)
    {
      insert(*it);
#if RMER_LOAD_STEP > 1
      for (unsigned i = 0; (i < RMER_LOAD_STEP) && (it != ntHashIterator::end()); ++i, ++it, ++offset) {
        if (offset == sequence.size() - m_kmerSize) { insert(*it); }
      }
#else
      ++it;
#endif
    }
  }
}

SpacedSeedsBloomFilter::SpacedSeedsBloomFilter(size_t filterSize, int kmerSize):
    BloomFilter(filterSize, SPACED_SEEDS_COUNT * SPACED_SEEDS_HASH_PER_SEED, kmerSize)
{
  const auto patterns = generateSpacedSeedsPatterns(SPACED_SEEDS_COUNT, kmerSize, SPACED_SEEDS_MISSES);
  for (const auto &pattern : patterns) {
    assert(pattern.size() == getKmerSize());
  }
  if (SPACED_SEEDS_QC) {
    QCSpacedSeedsPatterns(patterns);
  }
  spacedSeeds = stHashIterator::parseSeed(patterns);
}

size_t SpacedSeedsBloomFilter::getPop() const {
  size_t i, popBF = 0;
  #pragma omp parallel for reduction(+:popBF)
  for (i = 0; i < (m_size + 7) / 8; i++) {
    popBF = popBF + popCnt(m_filter[i]);
  }
  return popBF;
}

double SpacedSeedsBloomFilter::getFPR() const {
  const double occupancy = double(getPop()) / double(m_size);
  const double singleSeedFpr = std::pow(occupancy, SPACED_SEEDS_HASH_PER_SEED);
  const double totalFpr = 1 - std::pow(1 - singleSeedFpr, SPACED_SEEDS_COUNT);
  return totalFpr;
}

void SpacedSeedsBloomFilter::loadSequence(const Sequence& sequence) {
  if (sequence.size() >= m_kmerSize) {
#if RMER_LOAD_STEP > 1
    unsigned offset = 0;
#endif
    for (stHashIterator it(sequence, spacedSeeds, SPACED_SEEDS_COUNT, SPACED_SEEDS_HASH_PER_SEED, m_kmerSize);
      it != stHashIterator::end();)
    {
      insert(*it);
#if RMER_LOAD_STEP > 1
      for (unsigned i = 0; (i < RMER_LOAD_STEP) && (it != stHashIterator::end()); ++i, ++it, ++offset) {
        if (offset == sequence.size() - m_kmerSize) { insert(*it); }
      }
#else
      ++it;
#endif
    }
  }
}

static void loadReads(const std::vector<std::string>& readFilepaths, int r) {
  const size_t LOAD_PROGRESS_STEP = 100000;
  const size_t PARALLEL_IO_SIZE = 10;

  size_t kmersPerReadMode = 0;
  size_t kmersPerReadModeCount = 0;
  Histogram kmersPerReadHist;

  size_t threads_per_task = omp_get_max_threads() / PARALLEL_IO_SIZE;
  if (threads_per_task < 1) { threads_per_task = 1; }
  omp_set_nested(1);
  #pragma omp parallel
  #pragma omp single
  {
    int counter = 0;
    for (const auto path : readFilepaths)
    {
      #pragma omp task firstprivate(path)
      {
        assert(!path.empty());
        if (opt::verbose) {
          #pragma omp critical (cerr)
          std::cerr << "Loading reads from `" << path << "'...\n";
        }

        btllib::SeqReader reader(path);
        uint64_t readCount = 0, dropped = 0;
        #pragma omp parallel num_threads(threads_per_task)
        for (btllib::SeqReader::Record record; record = reader.read();) {
          if (int(record.seq.size()) != ReadBatch::current.size) { continue; }
          bool loaded = false;
          size_t qual_threshold_position = 0;
          for (int j = record.qual.size() - 1; j >= 0; j--) {
            if (record.qual[j] >= RMER_QUALITY_THRESHOLD) {
              qual_threshold_position = j;
              break;
            }
          }
          size_t substr_len = std::min(long(r + opt::threshold - 1),
            std::max(long(r + opt::threshold - 1), long(qual_threshold_position + 1)));
          std::string seq = record.seq.substr(0, substr_len);
          if (seq.size() >= g_vanillaBloom->getKmerSize()) {
            #pragma omp critical (kmersPerReadHist)
            kmersPerReadHist.insert(seq.size() - r + 1);

            g_vanillaBloom->loadSequence(seq);
            if (opt::error_correction) {
              g_spacedSeedsBloom->loadSequence(seq);
            }
            loaded = true;
          }
          if (opt::verbose)
          #pragma omp critical (cerr)
          {
            if (!loaded) { dropped++; }
            readCount++;
            if (readCount % LOAD_PROGRESS_STEP == 0)
              std::cerr << "\rLoaded " << readCount
                << " reads into Bloom filter.";
          }
        }
        if (opt::verbose) {
          std::cerr << std::endl;
        }
      }
      counter++;
      if (counter % PARALLEL_IO_SIZE == 0) {
        #pragma omp taskwait
      }
    }
  }

  for (const auto& entry : kmersPerReadHist) {
    if (entry.second > kmersPerReadModeCount) {
      kmersPerReadMode = entry.first;
      kmersPerReadModeCount = entry.second;
    }
  }
  if (opt::verbose) {
    std::cerr << "Kmers per read mode: " << kmersPerReadMode << '\n';
  }
}

VanillaBloomFilter* g_vanillaBloom = nullptr;
SpacedSeedsBloomFilter* g_spacedSeedsBloom = nullptr;

void buildFilters(const std::vector<std::string>& readFilepaths, const int r, const size_t bloomBytesTotal)
{
  assert(bloomBytesTotal > 0);
  try {
      if (opt::verbose) { std::cerr << "Building Bloom filter(s) for r value " << r << '\n'; }

      delete g_vanillaBloom;
      delete g_spacedSeedsBloom;
      
      size_t bloomBitsVanilla = size_t(bloomBytesTotal) * 8;
      size_t bloomBitsSpacedSeeds = 0;

      if (opt::error_correction) {
        double vanillaRatio = VANILLA_TO_SEEDS_MEM_RATIO * double(HASH_NUM) / (double(HASH_NUM) + double(SPACED_SEEDS_COUNT * SPACED_SEEDS_HASH_PER_SEED));
        bloomBitsVanilla = size_t(vanillaRatio * bloomBytesTotal) * 8;
        bloomBitsSpacedSeeds = bloomBytesTotal * 8 - bloomBitsVanilla;
      }
      

      if (opt::verbose > 1) {
        if (opt::error_correction) {
          std::cerr << "Total Bloom filter memory = " << bytesToSI(bloomBytesTotal) << '\n';
          std::cerr << "Vanilla Bloom filter memory = " << bytesToSI(bloomBitsVanilla / 8) << '\n';
          std::cerr << "Spaced seeds Bloom filter memory = " << bytesToSI(bloomBitsSpacedSeeds / 8) << '\n';
        } else {
          std::cerr << "Vanilla Bloom filter memory = " << bytesToSI(bloomBitsVanilla / 8) << '\n';
        }
      }

      g_vanillaBloom = new VanillaBloomFilter(bloomBitsVanilla, r);
      if (opt::error_correction) {
        g_spacedSeedsBloom = new SpacedSeedsBloomFilter(bloomBitsSpacedSeeds, r);
      }

      loadReads(readFilepaths, r);

      if (opt::verbose > 1) {
        const auto vanillaFPR = g_vanillaBloom->getFPR();

        std::cerr << "Vanilla Bloom filter (k = " << std::to_string(g_vanillaBloom->getKmerSize()) << setprecision(3)
          << ") occupancy = " << double(g_vanillaBloom->getPop()) / double(g_vanillaBloom->getFilterSize()) * 100.0 << "%"
          << ", FPR = " << vanillaFPR * 100.0 << "%" << endl;

        if (opt::error_correction) {
          std::cerr << "FPR for base substitution = " <<
            (1 - std::pow(1 - vanillaFPR, 3 * g_vanillaBloom->getKmerSize() / SPACED_SEEDS_COUNT * SPACED_SEEDS_SNP_FRACTION)) * 100.0 << "%" << std::endl;

          std::cerr << "Spaced seeds Bloom filter (k = " << std::to_string(g_spacedSeedsBloom->getKmerSize()) << setprecision(3)
            << ") occupancy = " << double(g_spacedSeedsBloom->getPop()) / double(g_spacedSeedsBloom->getFilterSize()) * 100.0 << "%"
            << ", FPR = " << g_spacedSeedsBloom->getFPR() * 100.0 << "%" << endl;
        }
      }
  } catch (const std::bad_alloc& e) {
    std::cerr << "Bloom filter allocation failed: " << e.what() << '\n';
    exit(EXIT_FAILURE);
  }
}
