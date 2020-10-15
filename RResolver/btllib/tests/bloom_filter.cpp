#include "../include/btllib/bloom_filter.hpp"

#include "helpers.hpp"

#include <cassert>
#include <cstdio>
#include <iostream>
#include <mutex>
#include <string>

int
main()
{
  std::cerr << "Testing BloomFilter" << std::endl;
  btllib::BloomFilter bf(1024 * 1024, 3);
  bf.insert({ 1, 10, 100 });
  bf.insert({ 100, 200, 300 });

  assert(bf.contains({ 1, 10, 100 }));
  assert(bf.contains({ 100, 200, 300 }));
  assert(!bf.contains({ 1, 20, 100 }));

  auto filename = get_random_name(64);
  bf.write(filename);

  btllib::BloomFilter bf2(filename);

  assert(bf2.contains({ 1, 10, 100 }));
  assert(bf2.contains({ 100, 200, 300 }));
  assert(!bf2.contains({ 1, 20, 100 }));

  std::remove(filename.c_str());

  std::string seq = "CACTATCGACGATCATTCGAGCATCAGCGACTG";
  std::string seq2 = "GTAGTACGATCAGCGACTATCGAGCTACGAGCA";
  assert(seq.size() == seq2.size());

  btllib::KmerBloomFilter kmer_bf(1024 * 1024, 4, seq.size() / 2);
  kmer_bf.insert(seq);
  assert(kmer_bf.contains(seq) == (seq.size() - seq.size() / 2 + 1));
  assert(kmer_bf.contains(seq2) <= 1);

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
  assert(false_positives < 10);

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
      assert(kmer_bf2.contains(seq));
    }
  }

  return 0;
}