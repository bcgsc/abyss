#include "btllib/bloom_filter.hpp"

#include "helpers.hpp"

#include <algorithm>
#include <cstdio>
#include <iostream>
#include <mutex>
#include <string>

int
main()
{
  std::cerr << "Testing BloomFilter" << std::endl;
  btllib::BloomFilter bf(1024 * 1024, 3, "ntHash");
  bf.insert({ 1, 10, 100 });
  bf.insert({ 100, 200, 300 });

  TEST_ASSERT(bf.contains({ 1, 10, 100 }));
  TEST_ASSERT(bf.contains({ 100, 200, 300 }));
  TEST_ASSERT(!bf.contains({ 1, 20, 100 }));

  auto filename = get_random_name(64);
  bf.save(filename);

  btllib::BloomFilter bf2(filename);
  TEST_ASSERT_EQ(bf2.get_hash_fn(), "ntHash");

  TEST_ASSERT(bf2.contains({ 1, 10, 100 }));
  TEST_ASSERT(bf2.contains({ 100, 200, 300 }));
  TEST_ASSERT(!bf2.contains({ 1, 20, 100 }));

  TEST_ASSERT(!bf2.contains_insert({ 9, 99, 999 }));
  TEST_ASSERT(bf2.contains_insert({ 9, 99, 999 }));

  std::remove(filename.c_str());

  std::string seq = "CACTATCGACGATCATTCGAGCATCAGCGACTG";
  std::string seq2 = "GTAGTACGATCAGCGACTATCGAGCTACGAGCA";
  TEST_ASSERT_EQ(seq.size(), seq2.size());

  std::cerr << "Testing KmerBloomFilter" << std::endl;
  btllib::KmerBloomFilter kmer_bf(1024 * 1024, 4, seq.size() / 2);
  kmer_bf.insert(seq);
  TEST_ASSERT_EQ(kmer_bf.contains(seq), (seq.size() - seq.size() / 2 + 1));
  TEST_ASSERT_LE(kmer_bf.contains(seq2), 1);

  std::cerr << "Testing SeedBloomFilter" << std::endl;
  std::string seed1 = "000001111111111111111111111111111";
  std::string seed2 = "111111111111111111111111111100000";
  std::string snp_seq1 = "AACTATCGACGATCATTCGAGCATCAGCGACTG";
  std::string snp_seq2 = "CACTATCGACGATCATTCGAGCATCAGCGACTA";
  TEST_ASSERT_EQ(seed1.size(), seed2.size());
  btllib::SeedBloomFilter seed_bf(1024 * 1024, seq.size(), { seed1, seed2 }, 4);
  seed_bf.insert(seq);
  auto hit_seeds = seed_bf.contains(seq);
  TEST_ASSERT(std::find(hit_seeds[0].begin(), hit_seeds[0].end(), 0) !=
              hit_seeds[0].end());
  TEST_ASSERT(std::find(hit_seeds[0].begin(), hit_seeds[0].end(), 1) !=
              hit_seeds[0].end());
  hit_seeds = seed_bf.contains(snp_seq1);
  TEST_ASSERT(std::find(hit_seeds[0].begin(), hit_seeds[0].end(), 0) !=
              hit_seeds[0].end());
  TEST_ASSERT(std::find(hit_seeds[0].begin(), hit_seeds[0].end(), 1) ==
              hit_seeds[0].end());
  hit_seeds = seed_bf.contains(snp_seq2);
  TEST_ASSERT(std::find(hit_seeds[0].begin(), hit_seeds[0].end(), 0) ==
              hit_seeds[0].end());
  TEST_ASSERT(std::find(hit_seeds[0].begin(), hit_seeds[0].end(), 1) !=
              hit_seeds[0].end());

  std::cerr << "Testing KmerBloomFilter with multiple threads" << std::endl;

  std::vector<std::string> present_seqs;
  std::vector<std::string> absent_seqs;
  for (size_t i = 0; i < 1000; i++) {
    present_seqs.push_back(get_random_seq(get_random(100, 200)));
    absent_seqs.push_back(get_random_seq(get_random(100, 200)));
  }
  std::vector<std::string> present_seqs2 = present_seqs;

  std::mutex seqs_lock;
  btllib::KmerBloomFilter kmer_bf2(50 * 1024 * 1024, 4, 100);
#pragma omp parallel shared(                                                   \
  present_seqs, present_seqs2, absent_seqs, seqs_lock, kmer_bf2)
  {
    while (true) {
      std::string seq;
      {
        std::unique_lock<std::mutex> lock(seqs_lock);
        if (present_seqs.empty()) {
          break;
        }
        seq = present_seqs.back();
        present_seqs.pop_back();
      }
      kmer_bf2.insert(seq);
    }
  }

  std::atomic<unsigned> false_positives(0);
#pragma omp parallel shared(present_seqs,                                      \
                            present_seqs2,                                     \
                            absent_seqs,                                       \
                            seqs_lock,                                         \
                            kmer_bf2,                                          \
                            false_positives)
  {
    while (true) {
      std::string seq;
      {
        std::unique_lock<std::mutex> lock(seqs_lock);
        if (absent_seqs.empty()) {
          break;
        }
        seq = absent_seqs.back();
        absent_seqs.pop_back();
      }
      false_positives += kmer_bf2.contains(seq);
    }
  }
  std::cerr << "False positives = " << false_positives << std::endl;
  TEST_ASSERT_LT(false_positives, 10);

#pragma omp parallel shared(                                                   \
  present_seqs, present_seqs2, absent_seqs, seqs_lock, kmer_bf2)
  {
    while (true) {
      std::string seq;
      {
        std::unique_lock<std::mutex> lock(seqs_lock);
        if (present_seqs2.empty()) {
          break;
        }
        seq = present_seqs2.back();
        present_seqs2.pop_back();
      }
      TEST_ASSERT(kmer_bf2.contains(seq));
    }
  }

  return 0;
}