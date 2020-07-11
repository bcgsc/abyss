#include "../include/btllib/bloom_filter.hpp"

#include "helpers.hpp"

#include <cassert>
#include <cstdio>
#include <iostream>
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

  std::cerr << "Testing CountingBloomFilter" << std::endl;
  btllib::CountingBloomFilter8 cbf(1024 * 1024, 3);

  cbf.insert({ 1, 10, 100 });
  cbf.insert({ 1, 10, 100 });
  cbf.insert({ 100, 200, 300 });

  assert(cbf.contains({ 1, 10, 100 }) == 2);
  assert(cbf.contains({ 100, 200, 300 }) == 1);
  assert(cbf.contains({ 1, 20, 100 }) == 0);

  auto filename2 = get_random_name(64);
  cbf.write(filename2);

  btllib::CountingBloomFilter8 cbf2(filename2);

  assert(cbf2.contains({ 1, 10, 100 }) == 2);
  assert(cbf2.contains({ 100, 200, 300 }) == 1);
  assert(cbf2.contains({ 1, 20, 100 }) == 0);

  std::remove(filename2.c_str());

  std::string seq = "CACTATCGACGATCATTCGAGCATCAGCGACTG";
  std::string seq2 = "GTAGTACGATCAGCGACTATCGAGCTACGAGCA";
  assert(seq.size() == seq2.size());

  btllib::KmerBloomFilter kmer_bf(seq.size() / 2, 1024 * 1024);
  kmer_bf.insert(seq);
  assert(kmer_bf.contains(seq) == (seq.size() - seq.size() / 2 + 1));
  assert(kmer_bf.contains(seq2) <= 1);

  return 0;
}