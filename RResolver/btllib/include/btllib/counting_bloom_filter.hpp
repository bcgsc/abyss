#ifndef BTLLIB_COUNTING_BLOOM_FILTER_HPP
#define BTLLIB_COUNTING_BLOOM_FILTER_HPP

#include "bloom_filter.hpp"
#include "nthash.hpp"
#include "status.hpp"

#include "../external/cpptoml.hpp"

#include <atomic>
#include <climits>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <string>
#include <vector>

namespace btllib {

static const char* const COUNTING_BLOOM_FILTER_MAGIC_HEADER =
  "BTLCountingBloomFilter_v5";
static const char* const KMER_COUNTING_BLOOM_FILTER_MAGIC_HEADER =
  "BTLKmerCountingBloomFilter_v5";

template<typename T>
class KmerCountingBloomFilter;

/**
 * Counting Bloom filter data structure. Provides CountingBloomFilter8,
 * CountingBloomFilter16, and CountingBloomFilter32 classes with corresponding
 * bit-size counters.
 */
template<typename T>
class CountingBloomFilter
{

public:
  /** Construct a dummy Kmer Bloom filter (e.g. as a default argument). */
  CountingBloomFilter() {}

  /**
   * Construct an empty Counting Bloom filter of given size.
   *
   * @param bytes Filter size in bytes.
   * @param hash_num Number of hash values per element.
   * @param hash_fn Name of the hash function used. Used for metadata. Optional.
   */
  CountingBloomFilter(size_t bytes,
                      unsigned hash_num,
                      std::string hash_fn = "");

  /**
   * Load a Counting Bloom filter from a file.
   *
   * @param path Filepath to load from.
   */
  explicit CountingBloomFilter(const std::string& path);

  ~CountingBloomFilter() { delete[] array; }

  CountingBloomFilter(const CountingBloomFilter&) = delete;
  CountingBloomFilter(CountingBloomFilter&&) = delete;

  CountingBloomFilter& operator=(const CountingBloomFilter&) = delete;
  CountingBloomFilter& operator=(CountingBloomFilter&&) = delete;

  /**
   * Insert an element's hash values.
   *
   * @param hashes Integer array of hash values. Array size should equal the
   * hash_num argument used when the Bloom filter was constructed.
   */
  void insert(const uint64_t* hashes);

  /**
   * Insert an element's hash values.
   *
   * @param hashes Integer vector of hash values.
   */
  void insert(const std::vector<uint64_t>& hashes) { insert(hashes.data()); }

  /**
   * Check for the presence of an element's hash values.
   *
   * @param hashes Integer array of hash values. Array size should equal the
   * hash_num argument used when the Bloom filter was constructed.
   *
   * @return The count of the queried element.
   */
  T contains(const uint64_t* hashes) const;

  /**
   * Check for the presence of an element's hash values.
   *
   * @param hashes Integer vector of hash values.
   *
   * @return The count of the queried element.
   */
  T contains(const std::vector<uint64_t>& hashes) const
  {
    return contains(hashes.data());
  }

  /**
   * Check for the presence of an element's hash values and insert if missing.
   *
   * @param hashes Integer array of hash values. Array size should equal the
   * hash_num argument used when the Bloom filter was constructed.
   *
   * @return The count of the queried element before insertion.
   */
  T contains_insert(const uint64_t* hashes);

  /**
   * Check for the presence of an element's hash values and insert if missing.
   *
   * @param hashes Integer vector of hash values.
   *
   * @return The count of the queried element before insertion.
   */
  T contains_insert(const std::vector<uint64_t>& hashes)
  {
    return contains_insert(hashes.data());
  }

  /** Get filter size in bytes. */
  size_t get_bytes() const { return bytes; }
  /** Get population count, i.e. the number of counters >0 in the filter. */
  uint64_t get_pop_cnt() const;
  /** Get the fraction of the filter occupied by >1 counters. */
  double get_occupancy() const;
  /** Get the number of hash values per element. */
  unsigned get_hash_num() const { return hash_num; }
  /** Get the query false positive rate. */
  double get_fpr() const;
  /** Get the name of the hash function used. */
  const std::string& get_hash_fn() const { return hash_fn; }

  /**
   * Save the Bloom filter to a file that can be loaded in the future.
   *
   * @param path Filepath to store filter at.
   */
  void save(const std::string& path);

private:
  friend class KmerCountingBloomFilter<T>;

  std::atomic<T>* array = nullptr;
  size_t bytes = 0;
  size_t array_size = 0;
  unsigned hash_num = 0;
  std::string hash_fn;
};

/**
 * Counting Bloom filter data structure that stores k-mers. Provides
 * KmerCountingBloomFilter8, KmerCountingBloomFilter16, and
 * KmerCountingBloomFilter32 classes with corresponding bit-size counters.
 */
template<typename T>
class KmerCountingBloomFilter
{

public:
  /** Construct a dummy Kmer Bloom filter (e.g. as a default argument). */
  KmerCountingBloomFilter() {}

  /**
   * Construct an empty Kmer Counting Bloom filter of given size.
   *
   * @param bytes Filter size in bytes.
   * @param hash_num Number of hash values per element.
   * @param k K-mer size.
   */
  KmerCountingBloomFilter(size_t bytes, unsigned hash_num, unsigned k);

  /**
   * Load a Kmer Counting Bloom filter from a file.
   *
   * @param path Filepath to load from.
   */
  explicit KmerCountingBloomFilter(const std::string& path);

  KmerCountingBloomFilter(const KmerCountingBloomFilter&) = delete;
  KmerCountingBloomFilter(KmerCountingBloomFilter&&) = delete;

  KmerCountingBloomFilter& operator=(const KmerCountingBloomFilter&) = delete;
  KmerCountingBloomFilter& operator=(KmerCountingBloomFilter&&) = delete;

  /**
   * Insert a sequence's k-mers into the filter.
   *
   * @param seq Sequence to k-merize.
   * @param seq_len Length of seq.
   */
  void insert(const char* seq, size_t seq_len);

  /**
   * Insert a sequence's k-mers into the filter.
   *
   * @param seq Sequence to k-merize.
   */
  void insert(const std::string& seq) { insert(seq.c_str(), seq.size()); }

  /**
   * Insert an element's hash values.
   *
   * @param hashes Integer array of hash values. Array size should equal the
   * hash_num argument used when the Bloom filter was constructed.
   */
  void insert(const uint64_t* hashes) { counting_bloom_filter.insert(hashes); }

  /**
   * Insert an element's hash values.
   *
   * @param hashes Integer vector of hash values.
   */
  void insert(const std::vector<uint64_t>& hashes)
  {
    counting_bloom_filter.insert(hashes.data());
  }

  /**
   * Query the presence of k-mers of a sequence.
   *
   * @param seq Sequence to k-merize.
   * @param seq_len Length of seq.
   *
   * @return The sum of counters of seq's k-mers found in the filter.
   */
  uint64_t contains(const char* seq, size_t seq_len) const;

  /**
   * Query the presence of k-mers of a sequence.
   *
   * @param seq Sequence to k-merize.
   *
   * @return The sum of counters of seq's k-mers found in the filter.
   */
  uint64_t contains(const std::string& seq) const
  {
    return contains(seq.c_str(), seq.size());
  }

  /**
   * Check for the presence of an element's hash values.
   *
   * @param hashes Integer array of hash values. Array size should equal the
   * hash_num argument used when the Bloom filter was constructed.
   *
   * @return The count of the queried element.
   */
  T contains(const uint64_t* hashes) const
  {
    return counting_bloom_filter.contains(hashes);
  }

  /**
   * Check for the presence of an element's hash values.
   *
   * @param hashes Integer vector of hash values.
   *
   * @return The count of the queried element.
   */
  T contains(const std::vector<uint64_t>& hashes) const
  {
    return counting_bloom_filter.contains(hashes);
  }

  /**
   * Check for the presence of sequence k-mers and insert if missing.
   *
   * @param hashes Integer array of hash values. Array size should equal the
   * hash_num argument used when the Bloom filter was constructed.
   *
   * @return The count of the queried element before insertion.
   */
  T contains_insert(const char* seq, size_t seq_len);

  /**
   * Check for the presence of sequence k-mers and insert if missing.
   *
   * @param hashes Integer vector of hash values.
   *
   * @return The count of the queried element before insertion.
   */
  T contains_insert(const std::string& seq)
  {
    return contains_insert(seq.c_str(), seq.size());
  }

  /**
   * Check for the presence of an element's hash values and insert if missing.
   *
   * @param hashes Integer array of hash values. Array size should equal the
   * hash_num argument used when the Bloom filter was constructed.
   *
   * @return The count of the queried element before insertion.
   */
  T contains_insert(const uint64_t* hashes)
  {
    return counting_bloom_filter.contains_insert(hashes);
  }

  /**
   * Check for the presence of an element's hash values and insert if missing.
   *
   * @param hashes Integer vector of hash values.
   *
   * @return The count of the queried element before insertion.
   */
  T contains_insert(const std::vector<uint64_t>& hashes)
  {
    return counting_bloom_filter.contains_insert(hashes.data());
  }

  /** Get filter size in bytes. */
  size_t get_bytes() const { return counting_bloom_filter.get_bytes(); }
  /** Get population count, i.e. the number of counters >0 in the filter. */
  uint64_t get_pop_cnt() const { return counting_bloom_filter.get_pop_cnt(); }
  /** Get the fraction of the filter occupied by >0 counters. */
  double get_occupancy() const { return counting_bloom_filter.get_occupancy(); }
  /** Get the number of hash values per element. */
  unsigned get_hash_num() const { return counting_bloom_filter.get_hash_num(); }
  /** Get the query false positive rate. */
  double get_fpr() const { return counting_bloom_filter.get_fpr(); }
  /** Get the k-mer size used. */
  unsigned get_k() const { return k; }
  /** Get the name of the hash function used. */
  const std::string& get_hash_fn() const
  {
    return counting_bloom_filter.get_hash_fn();
  }
  /** Get a reference to the underlying vanilla Counting Bloom filter. */
  CountingBloomFilter<T>& get_counting_bloom_filter()
  {
    return counting_bloom_filter;
  }

  /**
   * Save the Bloom filter to a file that can be loaded in the future.
   *
   * @param path Filepath to store filter at.
   */
  void save(const std::string& path);

private:
  CountingBloomFilter<T> counting_bloom_filter;
  unsigned k = 0;
};

using CountingBloomFilter8 = CountingBloomFilter<uint8_t>;
using CountingBloomFilter16 = CountingBloomFilter<uint16_t>;
using CountingBloomFilter32 = CountingBloomFilter<uint32_t>;

using KmerCountingBloomFilter8 = KmerCountingBloomFilter<uint8_t>;
using KmerCountingBloomFilter16 = KmerCountingBloomFilter<uint16_t>;
using KmerCountingBloomFilter32 = KmerCountingBloomFilter<uint32_t>;

template<typename T>
inline CountingBloomFilter<T>::CountingBloomFilter(size_t bytes,
                                                   unsigned hash_num,
                                                   std::string hash_fn)
  : bytes(
      size_t(std::ceil(double(bytes) / sizeof(uint64_t)) * sizeof(uint64_t)))
  , array_size(get_bytes() / sizeof(array[0]))
  , hash_num(hash_num)
  , hash_fn(std::move(hash_fn))
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
inline T
CountingBloomFilter<T>::contains_insert(const uint64_t* hashes)
{
  const auto prev_count = contains(hashes);
  insert(hashes);
  return prev_count;
}

template<typename T>
inline uint64_t
CountingBloomFilter<T>::get_pop_cnt() const
{
  uint64_t pop_cnt = 0;
#pragma omp parallel for default(none) reduction(+ : pop_cnt)
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
  if (table->contains("hash_fn")) {
    hash_fn = *(table->get_as<std::string>("hash_fn"));
  }
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
CountingBloomFilter<T>::save(const std::string& path)
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
  if (!hash_fn.empty()) {
    header->insert("hash_fn", hash_fn);
  }
  header->insert("counter_bits", size_t(sizeof(array[0]) * CHAR_BIT));
  root->insert(COUNTING_BLOOM_FILTER_MAGIC_HEADER, header);
  file << *root << "[HeaderEnd]\n";
  for (unsigned i = 0; i < PLACEHOLDER_NEWLINES; i++) {
    if (i == 1) {
      file << "  <binary data>";
    }
    file << '\n';
  }

  file.write((char*)array, array_size * sizeof(array[0]));
}

template<typename T>
inline KmerCountingBloomFilter<T>::KmerCountingBloomFilter(size_t bytes,
                                                           unsigned hash_num,
                                                           unsigned k)
  : counting_bloom_filter(bytes, hash_num, HASH_FN)
  , k(k)
{}

template<typename T>
inline void
KmerCountingBloomFilter<T>::insert(const char* seq, size_t seq_len)
{
  NtHash nthash(seq, seq_len, get_hash_num(), get_k());
  while (nthash.roll()) {
    counting_bloom_filter.insert(nthash.hashes());
  }
}

template<typename T>
inline uint64_t
KmerCountingBloomFilter<T>::contains(const char* seq, size_t seq_len) const
{
  uint64_t count = 0;
  NtHash nthash(seq, seq_len, get_hash_num(), get_k());
  while (nthash.roll()) {
    count += counting_bloom_filter.contains(nthash.hashes());
  }
  return count;
}

template<typename T>
inline T
KmerCountingBloomFilter<T>::contains_insert(const char* seq, size_t seq_len)
{
  const auto prev_count = contains(seq, seq_len);
  insert(seq, seq_len);
  return prev_count;
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
                "KmerCountingBloomFilter: Atomic primitives take extra memory. "
                "KmerCountingBloomFilter will "
                "have less than " +
                  std::to_string(get_bytes()) + " for bit array.");
  counting_bloom_filter.array_size =
    get_bytes() / sizeof(counting_bloom_filter.array[0]);
  counting_bloom_filter.hash_num =
    *table->get_as<decltype(counting_bloom_filter.hash_num)>("hash_num");
  const std::string loaded_hash_fn = *(table->get_as<std::string>("hash_fn"));
  check_error(
    loaded_hash_fn != HASH_FN,
    "KmerCountingBloomFilter: loaded hash function (" + loaded_hash_fn +
      ") is different from the one used by default (" + HASH_FN + ").");
  counting_bloom_filter.hash_fn = loaded_hash_fn;
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
KmerCountingBloomFilter<T>::save(const std::string& path)
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
  header->insert("hash_fn", get_hash_fn());
  header->insert("hash_num", get_hash_num());
  header->insert("counter_bits",
                 size_t(sizeof(counting_bloom_filter.array[0]) * CHAR_BIT));
  header->insert("k", k);
  root->insert(KMER_COUNTING_BLOOM_FILTER_MAGIC_HEADER, header);
  file << *root << "[HeaderEnd]\n";
  for (unsigned i = 0; i < PLACEHOLDER_NEWLINES; i++) {
    if (i == 1) {
      file << "  <binary data>";
    }
    file << '\n';
  }

  file.write((char*)counting_bloom_filter.array,
             counting_bloom_filter.array_size *
               sizeof(counting_bloom_filter.array[0]));
}

} // namespace btllib

#endif