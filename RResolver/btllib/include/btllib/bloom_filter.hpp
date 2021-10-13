#ifndef BTLLIB_BLOOM_FILTER_HPP
#define BTLLIB_BLOOM_FILTER_HPP

#include "nthash.hpp"
#include "status.hpp"

#include "../external/cpptoml.hpp"

#include <atomic>
#include <climits>
#include <cmath>
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

static const char* const BLOOM_FILTER_MAGIC_HEADER = "BTLBloomFilter_v6";
static const char* const KMER_BLOOM_FILTER_MAGIC_HEADER =
  "BTLKmerBloomFilter_v6";
static const char* const SEED_BLOOM_FILTER_MAGIC_HEADER =
  "BTLSeedBloomFilter_v6";
static const char* const HASH_FN = "ntHash_v1";

static const unsigned MAX_HASH_VALUES = 1024;
static const unsigned PLACEHOLDER_NEWLINES = 50;

inline unsigned
pop_cnt_byte(uint8_t x)
{
  return ((0x876543210 >>                                              // NOLINT
           (((0x4332322132212110 >> ((x & 0xF) << 2)) & 0xF) << 2)) >> // NOLINT
          ((0x4332322132212110 >> (((x & 0xF0) >> 2)) & 0xF) << 2)) &  // NOLINT
         0xf;                                                          // NOLINT
}

class BloomFilterInitializer
{

public:
  BloomFilterInitializer(const std::string& path,
                         const std::string& magic_string)
    : ifs(path)
    , table(parse_header(ifs, magic_string))
  {}

  /** Parse a Bloom filter file header. Useful for implementing Bloom filter
   * variants. */
  static std::shared_ptr<cpptoml::table> parse_header(
    std::ifstream& file,
    const std::string& magic_string);

  std::ifstream ifs;
  std::shared_ptr<cpptoml::table> table;

  BloomFilterInitializer(const BloomFilterInitializer&) = delete;
  BloomFilterInitializer(BloomFilterInitializer&&) = default;

  BloomFilterInitializer& operator=(const BloomFilterInitializer&) = delete;
  BloomFilterInitializer& operator=(BloomFilterInitializer&&) = default;
};

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

private:
  BloomFilter(BloomFilterInitializer bfi);

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

private:
  KmerBloomFilter(BloomFilterInitializer bfi);

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

private:
  SeedBloomFilter(BloomFilterInitializer bfi);

  std::vector<std::string> seeds;
  std::vector<SpacedSeed> parsed_seeds;
  KmerBloomFilter kmer_bloom_filter;
};

inline BloomFilter::BloomFilter(size_t bytes,
                                unsigned hash_num,
                                std::string hash_fn)
  : bytes(
      size_t(std::ceil(double(bytes) / sizeof(uint64_t)) * sizeof(uint64_t)))
  , array_size(get_bytes() / sizeof(array[0]))
  , array_bits(array_size * CHAR_BIT)
  , hash_num(hash_num)
  , hash_fn(std::move(hash_fn))
  , array(new std::atomic<uint8_t>[array_size])
{
  check_error(bytes == 0, "BloomFilter: memory budget must be >0!");
  check_error(hash_num == 0, "BloomFilter: number of hash values must be >0!");
  check_error(hash_num > MAX_HASH_VALUES,
              "BloomFilter: number of hash values cannot be over 1024!");
  check_warning(
    sizeof(uint8_t) != sizeof(std::atomic<uint8_t>),
    "Atomic primitives take extra memory. BloomFilter will have less than " +
      std::to_string(bytes) + " for bit array.");
  std::memset((void*)array.get(), 0, array_size * sizeof(array[0]));
}

inline void
BloomFilter::insert(const uint64_t* hashes)
{
  for (unsigned i = 0; i < hash_num; ++i) {
    const auto normalized = hashes[i] % array_bits;
    array[normalized / CHAR_BIT] |= BIT_MASKS[normalized % CHAR_BIT];
  }
}

inline bool
BloomFilter::contains(const uint64_t* hashes) const
{
  for (unsigned i = 0; i < hash_num; ++i) {
    const auto normalized = hashes[i] % array_bits;
    const auto mask = BIT_MASKS[normalized % CHAR_BIT];
    if (!bool(array[normalized / CHAR_BIT] & mask)) {
      return false;
    }
  }
  return true;
}

inline bool
BloomFilter::contains_insert(const uint64_t* hashes)
{
  uint8_t found = 1;
  for (unsigned i = 0; i < hash_num; ++i) {
    const auto normalized = hashes[i] % array_bits;
    const auto bitpos = normalized % CHAR_BIT;
    const auto mask = BIT_MASKS[bitpos];
    found &= ((array[normalized / CHAR_BIT].fetch_or(mask) >> bitpos) & 1);
  }
  return bool(found);
}

inline uint64_t
BloomFilter::get_pop_cnt() const
{
  uint64_t pop_cnt = 0;
#pragma omp parallel for default(none) reduction(+ : pop_cnt)
  for (size_t i = 0; i < array_size; ++i) {
    pop_cnt += pop_cnt_byte(array[i]);
  }
  return pop_cnt;
}

inline double
BloomFilter::get_occupancy() const
{
  return double(get_pop_cnt()) / double(array_bits);
}

inline double
BloomFilter::get_fpr() const
{
  return std::pow(get_occupancy(), double(hash_num));
}

inline std::shared_ptr<cpptoml::table>
BloomFilterInitializer::parse_header(std::ifstream& file,
                                     const std::string& magic_string)
{
  const std::string magic_with_brackets = std::string("[") + magic_string + "]";

  std::string line;
  std::getline(file, line);
  if (line != magic_with_brackets) {
    log_error(
      std::string("Magic string does not match (likely version mismatch)\n") +
      "File magic string:\t" + line + "\n" + "Loader magic string:\t" +
      magic_with_brackets);
    std::exit(EXIT_FAILURE); // NOLINT(concurrency-mt-unsafe)
  }

  /* Read bloom filter line by line until it sees "[HeaderEnd]"
  which is used to mark the end of the header section and
  assigns the header to a char array*/
  std::string toml_buffer(line + '\n');
  bool header_end_found = false;
  while (bool(std::getline(file, line))) {
    toml_buffer.append(line + '\n');
    if (line == "[HeaderEnd]") {
      header_end_found = true;
      break;
    }
  }
  if (!header_end_found) {
    log_error("Pre-built bloom filter does not have the correct header end.");
    std::exit(EXIT_FAILURE); // NOLINT(concurrency-mt-unsafe)
  }
  for (unsigned i = 0; i < PLACEHOLDER_NEWLINES; i++) {
    std::getline(file, line);
  }

  // Send the char array to a stringstream for the cpptoml parser to parse
  std::istringstream toml_stream(toml_buffer);
  cpptoml::parser toml_parser(toml_stream);
  const auto header_config = toml_parser.parse();

  // Obtain header values from toml parser and assign them to class members
  return header_config->get_table(magic_string);
}

inline BloomFilter::BloomFilter(const std::string& path)
  : BloomFilter::BloomFilter(
      BloomFilterInitializer(path, BLOOM_FILTER_MAGIC_HEADER))
{}

inline BloomFilter::BloomFilter(BloomFilterInitializer bfi)
  : bytes(*(bfi.table->get_as<decltype(bytes)>("bytes")))
  , array_size(bytes / sizeof(array[0]))
  , array_bits(array_size * CHAR_BIT)
  , hash_num(*(bfi.table->get_as<decltype(hash_num)>("hash_num")))
  , hash_fn(bfi.table->contains("hash_fn")
              ? *(bfi.table->get_as<decltype(hash_fn)>("hash_fn"))
              : "")
  , array(new std::atomic<uint8_t>[array_size])
{
  check_warning(
    sizeof(uint8_t) != sizeof(std::atomic<uint8_t>),
    "Atomic primitives take extra memory. BloomFilter will have less than " +
      std::to_string(bytes) + " for bit array.");
  bfi.ifs.read((char*)array.get(),
               std::streamsize(array_size * sizeof(array[0])));
}

inline void
BloomFilter::save(const std::string& path,
                  const cpptoml::table& table,
                  const char* data,
                  const size_t n)
{
  std::ofstream ofs(path.c_str(), std::ios::out | std::ios::binary);

  ofs << table << "[HeaderEnd]\n";
  for (unsigned i = 0; i < PLACEHOLDER_NEWLINES; i++) {
    if (i == 1) {
      ofs << "  <binary data>";
    }
    ofs << '\n';
  }

  ofs.write(data, std::streamsize(n));
}

inline void
BloomFilter::save(const std::string& path)
{
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
    header->insert("hash_fn", get_hash_fn());
  }
  root->insert(BLOOM_FILTER_MAGIC_HEADER, header);

  save(path, *root, (char*)array.get(), array_size * sizeof(array[0]));
}

inline KmerBloomFilter::KmerBloomFilter(size_t bytes,
                                        unsigned hash_num,
                                        unsigned k)
  : k(k)
  , bloom_filter(bytes, hash_num, HASH_FN)
{}

inline void
KmerBloomFilter::insert(const char* seq, size_t seq_len)
{
  NtHash nthash(seq, seq_len, get_hash_num(), get_k());
  while (nthash.roll()) {
    bloom_filter.insert(nthash.hashes());
  }
}

inline unsigned
KmerBloomFilter::contains(const char* seq, size_t seq_len) const
{
  unsigned count = 0;
  NtHash nthash(seq, seq_len, get_hash_num(), get_k());
  while (nthash.roll()) {
    if (bloom_filter.contains(nthash.hashes())) {
      count++;
    }
  }
  return count;
}

inline unsigned
KmerBloomFilter::contains_insert(const char* seq, size_t seq_len)
{
  unsigned count = 0;
  NtHash nthash(seq, seq_len, get_hash_num(), get_k());
  while (nthash.roll()) {
    if (bloom_filter.contains_insert(nthash.hashes())) {
      count++;
    }
  }
  return count;
}

inline KmerBloomFilter::KmerBloomFilter(const std::string& path)
  : KmerBloomFilter::KmerBloomFilter(
      BloomFilterInitializer(path, KMER_BLOOM_FILTER_MAGIC_HEADER))
{}

inline KmerBloomFilter::KmerBloomFilter(BloomFilterInitializer bfi)
  : k(*(bfi.table->get_as<decltype(k)>("k")))
  , bloom_filter(std::move(bfi))
{
  check_error(bloom_filter.hash_fn != HASH_FN,
              "KmerBloomFilter: loaded hash function (" + bloom_filter.hash_fn +
                ") is different from the one used by default (" + HASH_FN +
                ").");
}

inline void
KmerBloomFilter::save(const std::string& path)
{
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
  header->insert("hash_fn", get_hash_fn());
  header->insert("k", get_k());
  root->insert(KMER_BLOOM_FILTER_MAGIC_HEADER, header);

  BloomFilter::save(path,
                    *root,
                    (char*)bloom_filter.array.get(),
                    bloom_filter.array_size * sizeof(bloom_filter.array[0]));
}

inline SeedBloomFilter::SeedBloomFilter(size_t bytes,
                                        unsigned k,
                                        const std::vector<std::string>& seeds,
                                        unsigned hash_num_per_seed)
  : seeds(seeds)
  , parsed_seeds(parse_seeds(seeds))
  , kmer_bloom_filter(bytes, hash_num_per_seed, k)
{
  for (const auto& seed : seeds) {
    check_error(k != seed.size(),
                "SeedBloomFilter: passed k (" + std::to_string(k) +
                  ") not equal to passed spaced seed size (" +
                  std::to_string(seed.size()) + ")");
  }
}

inline void
SeedBloomFilter::insert(const char* seq, size_t seq_len)
{
  SeedNtHash nthash(
    seq, seq_len, parsed_seeds, get_hash_num_per_seed(), get_k());
  while (nthash.roll()) {
    for (size_t s = 0; s < seeds.size(); s++) {
      kmer_bloom_filter.bloom_filter.insert(nthash.hashes() +
                                            s * get_hash_num_per_seed());
    }
  }
}

inline std::vector<std::vector<unsigned>>
SeedBloomFilter::contains(const char* seq, size_t seq_len) const
{
  std::vector<std::vector<unsigned>> hit_seeds;
  SeedNtHash nthash(
    seq, seq_len, parsed_seeds, get_hash_num_per_seed(), get_k());
  while (nthash.roll()) {
    hit_seeds.emplace_back();
    for (size_t s = 0; s < seeds.size(); s++) {
      if (kmer_bloom_filter.bloom_filter.contains(
            nthash.hashes() + s * get_hash_num_per_seed())) {
        hit_seeds.back().push_back(s);
      }
    }
  }
  return hit_seeds;
}

inline std::vector<std::vector<unsigned>>
SeedBloomFilter::contains_insert(const char* seq, size_t seq_len)
{
  std::vector<std::vector<unsigned>> hit_seeds;
  SeedNtHash nthash(
    seq, seq_len, parsed_seeds, get_hash_num_per_seed(), get_k());
  while (nthash.roll()) {
    hit_seeds.emplace_back();
    for (size_t s = 0; s < seeds.size(); s++) {
      if (kmer_bloom_filter.bloom_filter.contains_insert(
            nthash.hashes() + s * get_hash_num_per_seed())) {
        hit_seeds.back().push_back(s);
      }
    }
  }
  return hit_seeds;
}

inline double
SeedBloomFilter::get_fpr() const
{
  const double single_seed_fpr =
    std::pow(get_occupancy(), get_hash_num_per_seed());
  return 1 - std::pow(1 - single_seed_fpr, seeds.size());
}

inline SeedBloomFilter::SeedBloomFilter(const std::string& path)
  : SeedBloomFilter::SeedBloomFilter(
      BloomFilterInitializer(path, SEED_BLOOM_FILTER_MAGIC_HEADER))
{}

inline SeedBloomFilter::SeedBloomFilter(BloomFilterInitializer bfi)
  : seeds(*(bfi.table->get_array_of<std::string>("seeds")))
  , parsed_seeds(parse_seeds(seeds))
  , kmer_bloom_filter(std::move(bfi))
{}

inline void
SeedBloomFilter::save(const std::string& path)
{
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
  header->insert("hash_fn", get_hash_fn());
  header->insert("k", get_k());
  auto seeds_array = cpptoml::make_array();
  for (const auto& seed : seeds) {
    seeds_array->push_back(seed);
  }
  header->insert("seeds", seeds_array);
  root->insert(SEED_BLOOM_FILTER_MAGIC_HEADER, header);

  BloomFilter::save(path,
                    *root,
                    (char*)kmer_bloom_filter.bloom_filter.array.get(),
                    kmer_bloom_filter.bloom_filter.array_size *
                      sizeof(kmer_bloom_filter.bloom_filter.array[0]));
}

} // namespace btllib

#endif
