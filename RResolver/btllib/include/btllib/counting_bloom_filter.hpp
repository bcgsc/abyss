#ifndef BTLLIB_COUNTING_BLOOM_FILTER_HPP
#define BTLLIB_COUNTING_BLOOM_FILTER_HPP

#include "bloom_filter.hpp"
#include "nthash.hpp"
#include "status.hpp"

#include "vendor/cpptoml.hpp"

#include <atomic>
#include <climits>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <string>
#include <vector>

namespace btllib {

static const char* const COUNTING_BLOOM_FILTER_MAGIC_HEADER =
  "BTLCountingBloomFilter_v2";
static const char* const KMER_COUNTING_BLOOM_FILTER_MAGIC_HEADER =
  "BTLKmerCountingBloomFilter_v2";

template<typename T>
class KmerCountingBloomFilter;

template<typename T>
class CountingBloomFilter
{

public:
  CountingBloomFilter() {}
  CountingBloomFilter(size_t bytes, unsigned hash_num);
  explicit CountingBloomFilter(const std::string& path);

  ~CountingBloomFilter() { delete[] array; }

  CountingBloomFilter(const CountingBloomFilter&) = delete;
  CountingBloomFilter(CountingBloomFilter&&) = delete;

  CountingBloomFilter& operator=(const CountingBloomFilter&) = delete;
  CountingBloomFilter& operator=(CountingBloomFilter&&) = delete;

  void insert(const uint64_t* hashes);
  void insert(const std::vector<uint64_t>& hashes) { insert(hashes.data()); }

  T contains(const uint64_t* hashes) const;
  T contains(const std::vector<uint64_t>& hashes) const
  {
    return contains(hashes.data());
  }

  size_t get_bytes() const { return bytes; }
  uint64_t get_pop_cnt() const;
  double get_occupancy() const;
  unsigned get_hash_num() const { return hash_num; }
  double get_fpr() const;

  /**
   * Write bloom filter data to a file
   * @param path output filepath
   */
  void write(const std::string& path);

private:
  friend class KmerCountingBloomFilter<T>;

  std::atomic<T>* array = nullptr;
  size_t bytes = 0;
  size_t array_size = 0;
  unsigned hash_num = 0;
};

template<typename T>
class KmerCountingBloomFilter
{

public:
  KmerCountingBloomFilter(size_t bytes, unsigned hash_num, unsigned k);
  explicit KmerCountingBloomFilter(const std::string& path);

  KmerCountingBloomFilter(const KmerCountingBloomFilter&) = delete;
  KmerCountingBloomFilter(KmerCountingBloomFilter&&) = delete;

  KmerCountingBloomFilter& operator=(const KmerCountingBloomFilter&) = delete;
  KmerCountingBloomFilter& operator=(KmerCountingBloomFilter&&) = delete;

  void insert(const char* seq, size_t seq_len);
  void insert(const std::string& seq) { insert(seq.c_str(), seq.size()); }

  uint64_t contains(const char* seq, size_t seq_len) const;
  uint64_t contains(const std::string& seq) const
  {
    return contains(seq.c_str(), seq.size());
  }

  T contains(const uint64_t* hashes) const
  {
    counting_bloom_filter.contains(hashes);
  }
  T contains(const std::vector<uint64_t>& hashes) const
  {
    counting_bloom_filter.contains(hashes);
  }

  size_t get_bytes() const { return counting_bloom_filter.get_bytes(); }
  uint64_t get_pop_cnt() const { return counting_bloom_filter.get_pop_cnt(); }
  double get_occupancy() const { return counting_bloom_filter.get_occupancy(); }
  unsigned get_hash_num() const { return counting_bloom_filter.get_hash_num(); }
  double get_fpr() const { return counting_bloom_filter.get_fpr(); }
  unsigned get_k() const { return k; }

  void write(const std::string& path);

private:
  CountingBloomFilter<T> counting_bloom_filter;
  unsigned k;
};

using CountingBloomFilter8 = CountingBloomFilter<uint8_t>;
using CountingBloomFilter16 = CountingBloomFilter<uint16_t>;
using CountingBloomFilter32 = CountingBloomFilter<uint32_t>;

using KmerCountingBloomFilter8 = KmerCountingBloomFilter<uint8_t>;
using KmerCountingBloomFilter16 = KmerCountingBloomFilter<uint16_t>;
using KmerCountingBloomFilter32 = KmerCountingBloomFilter<uint32_t>;

template<typename T>
inline CountingBloomFilter<T>::CountingBloomFilter(size_t bytes,
                                                   unsigned hash_num)
  : bytes(std::ceil(bytes / sizeof(uint64_t)) * sizeof(uint64_t))
  , array_size(get_bytes() / sizeof(array[0]))
  , hash_num(hash_num)
{
  check_warning(sizeof(uint8_t) != sizeof(std::atomic<uint8_t>),
                "Atomic primitives take extra memory. CountingBloomFilter will "
                "have less than " +
                  std::to_string(bytes) + " for bit array.");
  array = new std::atomic<T>[array_size];
  std::memset((void*)array, 0, array_size * sizeof(array[0]));
}

template<typename T>
inline void
CountingBloomFilter<T>::insert(const uint64_t* hashes)
{
  // Update flag to track if increment is done on at least one counter
  bool update_done = false;
  T new_val;
  T min_val = contains(hashes);
  while (!update_done) {
    // Simple check to deal with overflow
    new_val = min_val + 1;
    if (min_val > new_val) {
      return;
    }
    for (size_t i = 0; i < hash_num; ++i) {
      decltype(min_val) temp_min_val = min_val;
      if (array[hashes[i] % array_size].compare_exchange_strong(temp_min_val,
                                                                new_val)) {
        update_done = true;
      }
    }
    // Recalculate minval because if increment fails, it needs a new minval to
    // use and if it doesnt hava a new one, the while loop runs forever.
    if (!update_done) {
      min_val = contains(hashes);
    }
  }
}

template<typename T>
inline T
CountingBloomFilter<T>::contains(const uint64_t* hashes) const
{
  T min = array[hashes[0] % array_size];
  for (size_t i = 1; i < hash_num; ++i) {
    const size_t idx = hashes[i] % array_size;
    if (array[idx] < min) {
      min = array[idx];
    }
  }
  return min;
}

template<typename T>
inline uint64_t
CountingBloomFilter<T>::get_pop_cnt() const
{
  uint64_t pop_cnt = 0;
#pragma omp parallel for reduction(+ : pop_cnt)
  for (size_t i = 0; i < array_size; ++i) {
    if (array[i] > 0) {
      ++pop_cnt;
    }
  }
  return pop_cnt;
}

template<typename T>
inline double
CountingBloomFilter<T>::get_occupancy() const
{
  return double(get_pop_cnt()) / double(array_size);
}

template<typename T>
inline double
CountingBloomFilter<T>::get_fpr() const
{
  return std::pow(get_occupancy(), double(hash_num));
}

template<typename T>
inline CountingBloomFilter<T>::CountingBloomFilter(const std::string& path)
{
  std::ifstream file(path);

  auto table =
    BloomFilter::parse_header(file, COUNTING_BLOOM_FILTER_MAGIC_HEADER);
  bytes = *table->get_as<decltype(bytes)>("bytes");
  check_warning(sizeof(uint8_t) != sizeof(std::atomic<uint8_t>),
                "Atomic primitives take extra memory. CountingBloomFilter will "
                "have less than " +
                  std::to_string(bytes) + " for bit array.");
  array_size = bytes / sizeof(array[0]);
  hash_num = *table->get_as<decltype(hash_num)>("hash_num");
  check_error(
    sizeof(array[0]) * CHAR_BIT != *table->get_as<size_t>("counter_bits"),
    "CountingBloomFilter" + std::to_string(sizeof(array[0]) * CHAR_BIT) +
      " tried to load a file of CountingBloomFilter" +
      std::to_string(*table->get_as<size_t>("counter_bits")));

  array = new std::atomic<T>[array_size];
  file.read((char*)array, array_size * sizeof(array[0]));
}

template<typename T>
inline void
CountingBloomFilter<T>::write(const std::string& path)
{
  std::ofstream file(path.c_str(), std::ios::out | std::ios::binary);

  /* Initialize cpptoml root table
    Note: Tables and fields are unordered
    Ordering of table is maintained by directing the table
    to the output stream immediately after completion  */
  auto root = cpptoml::make_table();

  /* Initialize bloom filter section and insert fields
      and output to ostream */
  auto header = cpptoml::make_table();
  header->insert("bytes", get_bytes());
  header->insert("hash_num", get_hash_num());
  header->insert("counter_bits", size_t(sizeof(array[0]) * CHAR_BIT));
  root->insert(COUNTING_BLOOM_FILTER_MAGIC_HEADER, header);
  file << *root << "[HeaderEnd]\n";

  file.write((char*)array, array_size * sizeof(array[0]));
}

template<typename T>
inline KmerCountingBloomFilter<T>::KmerCountingBloomFilter(size_t bytes,
                                                           unsigned hash_num,
                                                           unsigned k)
  : counting_bloom_filter(bytes, hash_num)
  , k(k)
{}

template<typename T>
inline void
KmerCountingBloomFilter<T>::insert(const char* seq, size_t seq_len)
{
  NtHash nthash(seq, seq_len, get_k(), get_hash_num());
  while (nthash.roll()) {
    counting_bloom_filter.insert(nthash.hashes());
  }
}

template<typename T>
inline uint64_t
KmerCountingBloomFilter<T>::contains(const char* seq, size_t seq_len) const
{
  uint64_t count = 0;
  NtHash nthash(seq, seq_len, get_k(), get_hash_num());
  while (nthash.roll()) {
    count += counting_bloom_filter.contains(nthash.hashes());
  }
  return count;
}

template<typename T>
inline KmerCountingBloomFilter<T>::KmerCountingBloomFilter(
  const std::string& path)
{
  std::ifstream file(path);

  auto table =
    BloomFilter::parse_header(file, KMER_COUNTING_BLOOM_FILTER_MAGIC_HEADER);
  counting_bloom_filter.bytes =
    *table->get_as<decltype(counting_bloom_filter.bytes)>("bytes");
  check_warning(sizeof(uint8_t) != sizeof(std::atomic<uint8_t>),
                "Atomic primitives take extra memory. CountingBloomFilter will "
                "have less than " +
                  std::to_string(get_bytes()) + " for bit array.");
  counting_bloom_filter.array_size =
    get_bytes() / sizeof(counting_bloom_filter.array[0]);
  counting_bloom_filter.hash_num =
    *table->get_as<decltype(counting_bloom_filter.hash_num)>("hash_num");
  k = *table->get_as<decltype(k)>("k");
  check_error(sizeof(T) * CHAR_BIT != *table->get_as<size_t>("counter_bits"),
              "CountingBloomFilter" + std::to_string(sizeof(T) * CHAR_BIT) +
                " tried to load a file of CountingBloomFilter" +
                std::to_string(*table->get_as<size_t>("counter_bits")));

  counting_bloom_filter.array =
    new std::atomic<T>[counting_bloom_filter.array_size];
  file.read((char*)counting_bloom_filter.array,
            counting_bloom_filter.array_size *
              sizeof(counting_bloom_filter.array[0]));
}

template<typename T>
inline void
KmerCountingBloomFilter<T>::write(const std::string& path)
{
  std::ofstream file(path.c_str(), std::ios::out | std::ios::binary);

  /* Initialize cpptoml root table
    Note: Tables and fields are unordered
    Ordering of table is maintained by directing the table
    to the output stream immediately after completion  */
  auto root = cpptoml::make_table();

  /* Initialize bloom filter section and insert fields
      and output to ostream */
  auto header = cpptoml::make_table();
  header->insert("bytes", get_bytes());
  header->insert("hash_num", get_hash_num());
  header->insert("counter_bits",
                 size_t(sizeof(counting_bloom_filter.array[0]) * CHAR_BIT));
  header->insert("k", k);
  root->insert(KMER_COUNTING_BLOOM_FILTER_MAGIC_HEADER, header);
  file << *root << "[HeaderEnd]\n";

  file.write((char*)counting_bloom_filter.array,
             counting_bloom_filter.array_size *
               sizeof(counting_bloom_filter.array[0]));
}

} // namespace btllib

#endif
