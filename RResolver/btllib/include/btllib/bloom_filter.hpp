#ifndef BTLLIB_BLOOM_FILTER_HPP
#define BTLLIB_BLOOM_FILTER_HPP

#include "btllib/nthash.hpp"

#include "cpptoml.h"

#include <atomic>
#include <climits>
#include <cstdint>
#include <fstream>
#include <memory>
#include <string>
#include <vector>

namespace btllib {

static const uint8_t BIT_MASKS[CHAR_BIT] = {
  // NOLINT
  0x01, 0x02, 0x04, 0x08, // NOLINT
  0x10, 0x20, 0x40, 0x80  // NOLINT
};

static const char* const BLOOM_FILTER_SIGNATURE = "[BTLBloomFilter_v6]";
static const char* const KMER_BLOOM_FILTER_SIGNATURE =
  "[BTLKmerBloomFilter_v6]";
static const char* const SEED_BLOOM_FILTER_SIGNATURE =
  "[BTLSeedBloomFilter_v6]";
static const char* const HASH_FN = NTHASH_FN_NAME;

static const unsigned MAX_HASH_VALUES = 1024;
static const unsigned PLACEHOLDER_NEWLINES = 50;

/// @cond HIDDEN_SYMBOLS
class BloomFilterInitializer
{

public:
  BloomFilterInitializer(const std::string& path, const std::string& signature)
    : path(path)
    , ifs(path)
    , table(parse_header(signature))
  {
  }

  static bool check_file_signature(std::ifstream& ifs,
                                   const std::string& expected_signature,
                                   std::string& file_signature);

  std::string path;
  std::ifstream ifs;
  std::shared_ptr<cpptoml::table> table;

  BloomFilterInitializer(const BloomFilterInitializer&) = delete;
  BloomFilterInitializer(BloomFilterInitializer&&) = default;

  BloomFilterInitializer& operator=(const BloomFilterInitializer&) = delete;
  BloomFilterInitializer& operator=(BloomFilterInitializer&&) = default;

private:
  /** Parse a Bloom filter file header. Useful for implementing Bloom filter
   * variants. */
  std::shared_ptr<cpptoml::table> parse_header(const std::string& signature);
};
/// @endcond

class BloomFilter
{

public:
  /** Construct a dummy Bloom filter (e.g. as a default argument). */
  BloomFilter() {}

  /**
   * Construct an empty Bloom filter of given size.
   *
   * @param bytes Filter size in bytes.
   * @param hash_num Number of hash values per element.
   * @param hash_fn Name of the hash function used. Used for metadata. Optional.
   */
  BloomFilter(size_t bytes, unsigned hash_num, std::string hash_fn = "");

  /**
   * Load a Bloom filter from a file.
   *
   * @param path Filepath to load from.
   */
  explicit BloomFilter(const std::string& path);

  BloomFilter(const BloomFilter&) = delete;
  BloomFilter(BloomFilter&&) = delete;

  BloomFilter& operator=(const BloomFilter&) = delete;
  BloomFilter& operator=(BloomFilter&&) = delete;

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
   * @return True if present, false otherwise.
   */
  bool contains(const uint64_t* hashes) const;

  /**
   * Check for the presence of an element's hash values.
   *
   * @param hashes Integer vector of hash values.
   *
   * @return True if present, false otherwise.
   */
  bool contains(const std::vector<uint64_t>& hashes) const
  {
    return contains(hashes.data());
  }

  /**
   * Check for the presence of an element's hash values and insert if missing.
   *
   * @param hashes Integer array of hash values. Array size should equal the
   * hash_num argument used when the Bloom filter was constructed.
   *
   * @return True if present before insertion, false otherwise.
   */
  bool contains_insert(const uint64_t* hashes);

  /**
   * Check for the presence of an element's hash values and insert if missing.
   *
   * @param hashes Integer vector of hash values.
   *
   * @return True if present before insertion, false otherwise.
   */
  bool contains_insert(const std::vector<uint64_t>& hashes)
  {
    return contains_insert(hashes.data());
  }

  /** Get filter size in bytes. */
  size_t get_bytes() const { return bytes; }
  /** Get population count, i.e. the number of 1 bits in the filter. */
  uint64_t get_pop_cnt() const;
  /** Get the fraction of the filter occupied by 1 bits. */
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

  static void save(const std::string& path,
                   const cpptoml::table& table,
                   const char* data,
                   size_t n);

  /**
   * Check whether the file at the given path is a saved Bloom filter.
   *
   * @param path Filepath to check.
   */
  static bool is_bloom_file(const std::string& path)
  {
    return check_file_signature(path, BLOOM_FILTER_SIGNATURE);
  }

  static bool check_file_signature(const std::string& path,
                                   const std::string& signature);

private:
  BloomFilter(const std::shared_ptr<BloomFilterInitializer>& bfi);

  friend class KmerBloomFilter;
  friend class SeedBloomFilter;

  size_t bytes = 0;
  size_t array_size =
    0; // Should be equal to bytes, but not guaranteed by standard
  size_t array_bits = 0;
  unsigned hash_num = 0;
  std::string hash_fn;
  std::unique_ptr<std::atomic<uint8_t>[]> array;
};

/**
 * Bloom filter data structure stores k-mers.
 */
class KmerBloomFilter
{

public:
  /** Construct a dummy Kmer Bloom filter (e.g. as a default argument). */
  KmerBloomFilter() {}

  /**
   * Construct an empty Kmer Bloom filter of given size.
   *
   * @param bytes Filter size in bytes.
   * @param hash_num Number of hash values per element.
   * @param k K-mer size.
   */
  KmerBloomFilter(size_t bytes, unsigned hash_num, unsigned k);

  /**
   * Load a Kmer Bloom filter from a file.
   *
   * @param path Filepath to load from.
   */
  explicit KmerBloomFilter(const std::string& path);

  KmerBloomFilter(const KmerBloomFilter&) = delete;
  KmerBloomFilter(KmerBloomFilter&&) = delete;

  KmerBloomFilter& operator=(const KmerBloomFilter&) = delete;
  KmerBloomFilter& operator=(KmerBloomFilter&&) = delete;

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
  void insert(const uint64_t* hashes) { bloom_filter.insert(hashes); }

  /**
   * Insert an element's hash values.
   *
   * @param hashes Integer vector of hash values.
   */
  void insert(const std::vector<uint64_t>& hashes)
  {
    bloom_filter.insert(hashes);
  }

  /**
   * Query the presence of k-mers of a sequence.
   *
   * @param seq Sequence to k-merize.
   * @param seq_len Length of seq.
   *
   * @return The number of seq's k-mers found in the filter.
   */
  unsigned contains(const char* seq, size_t seq_len) const;

  /**
   * Query the presence of k-mers of a sequence.
   *
   * @param seq Sequence to k-merize.
   *
   * @return The number of seq's k-mers found in the filter.
   */
  unsigned contains(const std::string& seq) const
  {
    return contains(seq.c_str(), seq.size());
  }

  /**
   * Check for the presence of an element's hash values.
   *
   * @param hashes Integer array of hash values. Array size should equal the
   * hash_num argument used when the Bloom filter was constructed.
   */
  bool contains(const uint64_t* hashes) const
  {
    return bloom_filter.contains(hashes);
  }

  /**
   * Check for the presence of an element's hash values.
   *
   * @param hashes Integer vector of hash values.
   */
  bool contains(const std::vector<uint64_t>& hashes) const
  {
    return bloom_filter.contains(hashes);
  }

  /**
   * Query the presence of k-mers of a sequence and insert if missing.
   *
   * @param seq Sequence to k-merize.
   * @param seq_len Length of seq.
   *
   * @return The number of seq's k-mers found in the filter before insertion.
   */
  unsigned contains_insert(const char* seq, size_t seq_len);

  /**
   * Query the presence of k-mers of a sequence and insert if missing.
   *
   * @param seq Sequence to k-merize.
   *
   * @return The number of seq's k-mers found in the filter before insertion.
   */
  unsigned contains_insert(const std::string& seq)
  {
    return contains_insert(seq.c_str(), seq.size());
  }

  /**
   * Check for the presence of an element's hash values and insert if missing.
   *
   * @param hashes Integer array of hash values. Array size should equal the
   * hash_num argument used when the Bloom filter was constructed.
   *
   * @return True if present before insertion, false otherwise.
   */
  bool contains_insert(const uint64_t* hashes)
  {
    return bloom_filter.contains_insert(hashes);
  }

  /**
   * Check for the presence of an element's hash values and insert if missing.
   *
   * @param hashes Integer vector of hash values.
   *
   * @return True if present before insertion, false otherwise.
   */
  bool contains_insert(const std::vector<uint64_t>& hashes)
  {
    return bloom_filter.contains_insert(hashes);
  }

  /** Get filter size in bytes. */
  size_t get_bytes() const { return bloom_filter.get_bytes(); }
  /** Get population count, i.e. the number of 1 bits in the filter. */
  uint64_t get_pop_cnt() const { return bloom_filter.get_pop_cnt(); }
  /** Get the fraction of the filter occupied by 1 bits. */
  double get_occupancy() const { return bloom_filter.get_occupancy(); }
  /** Get the number of hash values per element. */
  unsigned get_hash_num() const { return bloom_filter.get_hash_num(); }
  /** Get the query false positive rate. */
  double get_fpr() const { return bloom_filter.get_fpr(); }
  /** Get the k-mer size used. */
  unsigned get_k() const { return k; }
  /** Get the name of the hash function used. */
  const std::string& get_hash_fn() const { return bloom_filter.get_hash_fn(); }
  /** Get a reference to the underlying vanilla Bloom filter. */
  BloomFilter& get_bloom_filter() { return bloom_filter; }

  /**
   * Save the Bloom filter to a file that can be loaded in the future.
   *
   * @param path Filepath to store filter at.
   */
  void save(const std::string& path);

  /**
   * Check whether the file at the given path is a saved Kmer Bloom filter.
   *
   * @param path Filepath to check.
   */
  static bool is_bloom_file(const std::string& path)
  {
    return btllib::BloomFilter::check_file_signature(
      path, KMER_BLOOM_FILTER_SIGNATURE);
  }

private:
  KmerBloomFilter(const std::shared_ptr<BloomFilterInitializer>& bfi);

  friend class SeedBloomFilter;

  unsigned k = 0;
  BloomFilter bloom_filter;
};

/**
 * Bloom filter data structure that stores spaced seed k-mers.
 */
class SeedBloomFilter
{

public:
  /** Construct a dummy Seed Bloom filter (e.g. as a default argument). */
  SeedBloomFilter() {}

  /**
   * Construct an empty Seed Bloom filter of given size.
   *
   * @param bytes Filter size in bytes.
   * @param k K-mer size.
   * @param seeds A vector of spaced seeds in string format. 0s indicate ignored
   * and 1s indicate relevant bases.
   */
  SeedBloomFilter(size_t bytes,
                  unsigned k,
                  const std::vector<std::string>& seeds,
                  unsigned hash_num_per_seed);

  /**
   * Load a Seed Bloom filter from a file.
   *
   * @param path Filepath to load from.
   */
  explicit SeedBloomFilter(const std::string& path);

  SeedBloomFilter(const SeedBloomFilter&) = delete;
  SeedBloomFilter(SeedBloomFilter&&) = delete;

  SeedBloomFilter& operator=(const SeedBloomFilter&) = delete;
  SeedBloomFilter& operator=(SeedBloomFilter&&) = delete;

  /**
   * Insert a sequence's spaced seed k-mers into the filter.
   *
   * @param seq Sequence to k-merize.
   * @param seq_len Length of seq.
   */
  void insert(const char* seq, size_t seq_len);

  /**
   * Insert a sequence's spaced seed k-mers into the filter.
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
  void insert(const uint64_t* hashes) { kmer_bloom_filter.insert(hashes); }

  /**
   * Insert an element's hash values.
   *
   * @param hashes Integer vector of hash values.
   */
  void insert(const std::vector<uint64_t>& hashes)
  {
    kmer_bloom_filter.insert(hashes);
  }

  /**
   * Query the presence of spaced seed k-mers of a sequence.
   *
   * @param seq Sequence to k-merize.
   * @param seq_len Length of seq.
   *
   * @return A vector indicating which seeds had a hit for every k-mer. The
   * indices of the outer vector are indices of seq k-mers. The indices of inner
   * vector are indices of spaced seeds hit for that k-mer.
   */
  std::vector<std::vector<unsigned>> contains(const char* seq,
                                              size_t seq_len) const;

  /**
   * Query the presence of spaced seed k-mers of a sequence.
   *
   * @param seq Sequence to k-merize.
   *
   * @return A vector indicating which seeds had a hit for every k-mer. The
   * indices of the outer vector are indices of seq k-mers. The indices of inner
   * vector are indices of spaced seeds hit for that k-mer.
   */
  std::vector<std::vector<unsigned>> contains(const std::string& seq) const
  {
    return contains(seq.c_str(), seq.size());
  }

  /**
   * Check for the presence of an element's hash values. A single spaced seed is
   * an element here.
   *
   * @param hashes Integer array of hash values. Array size should equal the
   * hash_num_per_seed argument used when the Bloom filter was constructed.
   */
  bool contains(const uint64_t* hashes) const
  {
    return kmer_bloom_filter.contains(hashes);
  }

  /**
   * Check for the presence of an element's hash values. A single spaced seed is
   * an element here.
   *
   * @param hashes Integer vector of hash values.
   */
  bool contains(const std::vector<uint64_t>& hashes) const
  {
    return kmer_bloom_filter.contains(hashes);
  }

  /**
   * Query the presence of spaced seed k-mers of a sequence and insert if
   * missing.
   *
   * @param seq Sequence to k-merize.
   * @param seq_len Length of seq.
   *
   * @return A vector indicating which seeds had a hit for every k-mer before
   * insertion. The indices of the outer vector are indices of seq k-mers. The
   * indices of inner vector are indices of spaced seeds hit for that k-mer.
   */
  std::vector<std::vector<unsigned>> contains_insert(const char* seq,
                                                     size_t seq_len);

  /**
   * Query the presence of spaced seed k-mers of a sequence and insert if
   * missing.
   *
   * @param seq Sequence to k-merize.
   *
   * @return A vector indicating which seeds had a hit for every k-mer before
   * insertion. The indices of the outer vector are indices of seq k-mers. The
   * indices of inner vector are indices of spaced seeds hit for that k-mer.
   */
  std::vector<std::vector<unsigned>> contains_insert(const std::string& seq)
  {
    return contains_insert(seq.c_str(), seq.size());
  }

  /**
   * Check for the presence of an element's hash values and insert if missing. A
   * single spaced seed is an element here.
   *
   * @param hashes Integer array of hash values. Array size should equal the
   * hash_num_per_seed argument used when the Bloom filter was constructed.
   *
   * @return True if present before insertion, false otherwise.
   */
  bool contains_insert(const uint64_t* hashes)
  {
    return kmer_bloom_filter.contains_insert(hashes);
  }

  /**
   * Check for the presence of an element's hash values and insert if missing. A
   * single spaced seed is an element here.
   *
   * @param hashes Integer vector of hash values.
   *
   * @return True if present before insertion, false otherwise.
   */
  bool contains_insert(const std::vector<uint64_t>& hashes)
  {
    return kmer_bloom_filter.contains_insert(hashes);
  }

  /** Get filter size in bytes. */
  size_t get_bytes() const { return kmer_bloom_filter.get_bytes(); }
  /** Get population count, i.e. the number of 1 bits in the filter. */
  uint64_t get_pop_cnt() const { return kmer_bloom_filter.get_pop_cnt(); }
  /** Get the fraction of the filter occupied by 1 bits. */
  double get_occupancy() const { return kmer_bloom_filter.get_occupancy(); }
  /** Get the number of hash values per k-mer, i.e. total number of hash values
   * for all seeds. */
  unsigned get_total_hash_num() const
  {
    return get_hash_num_per_seed() * get_seeds().size();
  }
  /** Get the false positive rate of at least one seed falsely reporting a hit
   * per k-mer. */
  double get_fpr() const;
  /** Get the k-mer size used. */
  unsigned get_k() const { return kmer_bloom_filter.get_k(); }
  /** Get the seeds used in string format. */
  const std::vector<std::string>& get_seeds() const { return seeds; }
  /** Get the seeds used in parsed format. Parsed format is a vector of indices
   * of 0s in the seed. */
  const std::vector<SpacedSeed>& get_parsed_seeds() const
  {
    return parsed_seeds;
  }
  /** Get the number of hash values per element, i.e. seed. */
  unsigned get_hash_num_per_seed() const
  {
    return kmer_bloom_filter.get_hash_num();
  }
  /** Get the number of hash values per element, i.e. seed. */
  unsigned get_hash_num() const { return get_hash_num_per_seed(); }
  /** Get the name of the hash function used. */
  const std::string& get_hash_fn() const
  {
    return kmer_bloom_filter.get_hash_fn();
  }
  /** Get a reference to the underlying Kmer Bloom filter. */
  KmerBloomFilter& get_kmer_bloom_filter() { return kmer_bloom_filter; }

  /**
   * Save the Bloom filter to a file that can be loaded in the future.
   *
   * @param path Filepath to store filter at.
   */
  void save(const std::string& path);

  /**
   * Check whether the file at the given path is a saved Seed Bloom filter.
   *
   * @param path Filepath to check.
   */
  static bool is_bloom_file(const std::string& path)
  {
    return btllib::BloomFilter::check_file_signature(
      path, SEED_BLOOM_FILTER_SIGNATURE);
  }

private:
  SeedBloomFilter(const std::shared_ptr<BloomFilterInitializer>& bfi);

  std::vector<std::string> seeds;
  std::vector<SpacedSeed> parsed_seeds;
  KmerBloomFilter kmer_bloom_filter;
};

} // namespace btllib

#endif
