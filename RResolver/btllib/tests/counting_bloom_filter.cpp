#include "btllib/counting_bloom_filter.hpp"

#include "helpers.hpp"

#include <cassert>
#include <cstdio>
#include <iostream>
#include <mutex>
#include <string>

int
main()
{
  std::cerr << "Testing CountingBloomFilter" << std::endl;
  btllib::CountingBloomFilter8 cbf(1024 * 1024, 3);

  cbf.insert({ 1, 10, 100 });
  cbf.insert({ 1, 10, 100 });
  cbf.insert({ 100, 200, 300 });

  assert(cbf.contains({ 1, 10, 100 }) == 2);
  assert(cbf.contains({ 100, 200, 300 }) == 1);
  assert(cbf.contains({ 1, 20, 100 }) == 0);

  auto filename = get_random_name(64);
  cbf.write(filename);

  btllib::CountingBloomFilter8 cbf2(filename);

  assert(cbf2.contains({ 1, 10, 100 }) == 2);
  assert(cbf2.contains({ 100, 200, 300 }) == 1);
  assert(cbf2.contains({ 1, 20, 100 }) == 0);

  std::remove(filename.c_str());

  std::string seq = "CACTATCGACGATCATTCGAGCATCAGCGACTG";
  std::string seq2 = "GTAGTACGATCAGCGACTATCGAGCTACGAGCA";
  assert(seq.size() == seq2.size());

  std::cerr << "Testing KmerCountingBloomFilter" << std::endl;
  btllib::KmerCountingBloomFilter8 kmer_counting_bf(
    1024 * 1024, 4, seq.size() / 2);
  kmer_counting_bf.insert(seq);
  assert(kmer_counting_bf.contains(seq) == (seq.size() - seq.size() / 2 + 1));
  assert(kmer_counting_bf.contains(seq2) <= 1);

  std::cerr << "Testing KmerCountingBloomfilter with multiple threads"
            << std::endl;

  std::vector<std::string> present_seqs;
  std::vector<unsigned> inserts;
  std::vector<std::string> absent_seqs;
  for (size_t i = 0; i < 100; i++) {
    present_seqs.push_back(get_random_seq(100));
    inserts.push_back(get_random(1, 10));
    absent_seqs.push_back(get_random_seq(100));
  }
  std::vector<std::string> present_seqs2 = present_seqs;

  std::mutex seqs_lock;
  btllib::KmerCountingBloomFilter8 kmer_counting_bf2(100 * 1024 * 1024, 4, 100);
#pragma omp parallel shared(present_seqs,                                      \
                            present_seqs2,                                     \
                            inserts,                                           \
                            absent_seqs,                                       \
                            seqs_lock,                                         \
                            kmer_counting_bf2)
  {
    while (true) {
      std::string seq;
      unsigned insert_count;
      {
        std::unique_lock<std::mutex> lock(seqs_lock);
        if (present_seqs.empty()) {
          break;
        }
        seq = present_seqs.back();
        insert_count = inserts.back();
        present_seqs.pop_back();
        inserts.pop_back();
      }
      for (size_t i = 0; i < insert_count; i++) {
        kmer_counting_bf2.insert(seq);
      }
    }
  }

  std::atomic<unsigned> false_positives(0);
#pragma omp parallel shared(present_seqs,                                      \
                            present_seqs2,                                     \
                            inserts,                                           \
                            absent_seqs,                                       \
                            seqs_lock,                                         \
                            kmer_counting_bf2,                                 \
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
      false_positives += kmer_counting_bf2.contains(seq);
    }
  }
  std::cerr << "False positives = " << false_positives << std::endl;
  assert(false_positives < 10);

  std::atomic<int> more_than_1(0);
#pragma omp parallel shared(present_seqs,                                      \
                            present_seqs2,                                     \
                            inserts,                                           \
                            absent_seqs,                                       \
                            seqs_lock,                                         \
                            kmer_counting_bf2,                                 \
                            more_than_1)
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
      if (kmer_counting_bf2.contains(seq) > 1) {
        more_than_1++;
      }
      assert(kmer_counting_bf2.contains(seq));
    }
  }
  std::cerr << "Seqs with more than 1 presence = " << more_than_1 << std::endl;
  assert(more_than_1 > 5);

  return 0;
}