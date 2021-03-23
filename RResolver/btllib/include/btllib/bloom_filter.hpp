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

static const unsigned char BIT_MASKS[CHAR_BIT] = {
  // NOLINT
  0x01, 0x02, 0x04, 0x08, // NOLINT
  0x10, 0x20, 0x40, 0x80  // NOLINT
};

static const char* const BLOOM_FILTER_MAGIC_HEADER = "BTLBloomFilter_v2";
static const char* const KMER_BLOOM_FILTER_MAGIC_HEADER =
  "BTLKmerBloomFilter_v2";
static const char* const SEED_BLOOM_FILTER_MAGIC_HEADER =
  "BTLSeedBloomFilter_v2";

static const unsigned MAX_HASH_VALUES = 1024;

inline unsigned
pop_cnt_byte(uint8_t x)
{
  return ((0x876543210 >>                                              // NOLINT
           (((0x4332322132212110 >> ((x & 0xF) << 2)) & 0xF) << 2)) >> // NOLINT
          ((0x4332322132212110 >> (((x & 0xF0) >> 2)) & 0xF) << 2)) &  // NOLINT
         0xf;                                                          // NOLINT
}

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
   */
  BloomFilter(size_t bytes, unsigned hash_num);

  /**
   * Load a Bloom filter from a file.
   *
   * @param path Filepath to load from.
   */
  explicit BloomFilter(const std::string& path);

  ~BloomFilter() { delete[] array; }

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

  /**
   * Write the Bloom filter to a file that can be loaded in the future.
   *
   * @param path Filepath to store filter at.
   */
  void write(const std::string& path);

  /** Parse a Bloom filter file header. Useful for implementing Bloom filter
   * variants. */
  static std::shared_ptr<cpptoml::table> parse_header(
    std::ifstream& file,
    const std::string& magic_string);

private:
  friend class KmerBloomFilter;
  friend class SeedBloomFilter;

  std::atomic<uint8_t>* array = nullptr;
  size_t bytes = 0;
  size_t array_size =
    0; // Should be equal to bytes, but not guaranteed by standard
  size_t array_bits = 0;
  unsigned hash_num = 0;
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
  /** Get a reference to the underlying vanilla Bloom filter. */
  BloomFilter& get_bloom_filter() { return bloom_filter; }

  /**
   * Write the Bloom filter to a file that can be loaded in the future.
   *
   * @param path Filepath to store filter at.
   */
  void write(const std::string& path);

private:
  friend class SeedBloomFilter;

  BloomFilter bloom_filter;
  unsigned k = 0;
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

  /** Get filter size in bytes. */
  size_t get_bytes() const { return kmer_bloom_filter.get_bytes(); }
  /** Get population count, i.e. the number of 1 bits in the filter. */
  uint64_t get_pop_cnt() const { return kmer_bloom_filter.get_pop_cnt(); }
  /** Get the fraction of the filter occupied by 1 bits. */
  double get_occupancy() const { return kmer_bloom_filter.get_occupancy(); }
  /** Get the number of hash values per k-mer, i.e. total number of hash values
   * for all seeds. */
  unsigned get_hash_num() const
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
  /** Get a reference to the underlying Kmer Bloom filter. */
  KmerBloomFilter& get_kmer_bloom_filter() { return kmer_bloom_filter; }

  /**
   * Write the Bloom filter to a file that can be loaded in the future.
   *
   * @param path Filepath to store filter at.
   */
  void write(const std::string& path);

private:
  KmerBloomFilter kmer_bloom_filter;
  std::vector<std::string> seeds;
  std::vector<SpacedSeed> parsed_seeds;
};

inline BloomFilter::BloomFilter(size_t bytes, unsigned hash_num)
  : bytes(std::ceil(bytes / sizeof(uint64_t)) * sizeof(uint64_t))
  , array_size(get_bytes() / sizeof(array[0]))
  , array_bits(array_size * CHAR_BIT)
  , hash_num(hash_num)
{
  check_error(bytes == 0, "BloomFilter: memory budget must be >0!");
  check_error(hash_num == 0, "BloomFilter: number of hash values must be >0!");
  check_error(hash_num > MAX_HASH_VALUES,
              "BloomFilter: number of hash values cannot be over 1024!");
  check_warning(
    sizeof(uint8_t) != sizeof(std::atomic<uint8_t>),
    "Atomic primitives take extra memory. BloomFilter will have less than " +
      std::to_string(bytes) + " for bit array.");
  array = new std::atomic<uint8_t>[array_size];
  std::memset((void*)array, 0, array_size * sizeof(array[0]));
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
BloomFilter::parse_header(std::ifstream& file, const std::string& magic_string)
{
  const std::string magic_with_brackets = std::string("[") + magic_string + "]";

  std::string line;
  std::getline(file, line);
  if (line != magic_with_brackets) {
    log_error(
      std::string("Magic string does not match (likely version mismatch)\n") +
      "File magic string:\t" + line + "\n" + "Loader magic string:\t" +
      magic_with_brackets);
    std::exit(EXIT_FAILURE);
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
    std::exit(EXIT_FAILURE);
  }

  // Send the char array to a stringstream for the cpptoml parser to parse
  std::istringstream toml_stream(toml_buffer);
  cpptoml::parser toml_parser(toml_stream);
  const auto header_config = toml_parser.parse();

  // Obtain header values from toml parser and assign them to class members
  return header_config->get_table(magic_string);
}

inline BloomFilter::BloomFilter(const std::string& path)
{
  std::ifstream file(path);

  auto table = parse_header(file, BLOOM_FILTER_MAGIC_HEADER);
  bytes = *(table->get_as<decltype(bytes)>("bytes"));
  check_warning(
    sizeof(uint8_t) != sizeof(std::atomic<uint8_t>),
    "Atomic primitives take extra memory. BloomFilter will have less than " +
      std::to_string(bytes) + " for bit array.");
  array_size = bytes / sizeof(std::atomic<uint8_t>);
  array_bits = array_size * CHAR_BIT;
  hash_num = *(table->get_as<decltype(hash_num)>("hash_num"));

  array = new std::atomic<uint8_t>[array_size];
  file.read((char*)array, array_size * sizeof(array[0]));
}

inline void
BloomFilter::write(const std::string& path)
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
  root->insert(BLOOM_FILTER_MAGIC_HEADER, header);
  file << *root << "[HeaderEnd]\n";

  file.write((char*)array, array_size * sizeof(array[0]));
}

inline KmerBloomFilter::KmerBloomFilter(size_t bytes,
                                        unsigned hash_num,
                                        unsigned k)
  : bloom_filter(bytes, hash_num)
  , k(k)
{}

inline void
KmerBloomFilter::insert(const char* seq, size_t seq_len)
{
  NtHash nthash(seq, seq_len, get_k(), get_hash_num());
  while (nthash.roll()) {
    bloom_filter.insert(nthash.hashes());
  }
}

inline unsigned
KmerBloomFilter::contains(const char* seq, size_t seq_len) const
{
  unsigned count = 0;
  NtHash nthash(seq, seq_len, get_k(), get_hash_num());
  while (nthash.roll()) {
    if (bloom_filter.contains(nthash.hashes())) {
      count++;
    }
  }
  return count;
}

inline KmerBloomFilter::KmerBloomFilter(const std::string& path)
{
  std::ifstream file(path);

  auto table = bloom_filter.parse_header(file, KMER_BLOOM_FILTER_MAGIC_HEADER);
  bloom_filter.bytes = *(table->get_as<decltype(bloom_filter.bytes)>("bytes"));
  check_warning(
    sizeof(uint8_t) != sizeof(std::atomic<uint8_t>),
    "Atomic primitives take extra memory. BloomFilter will have less than " +
      std::to_string(get_bytes()) + " for bit array.");
  bloom_filter.array_size = get_bytes() / sizeof(bloom_filter.array[0]);
  bloom_filter.array_bits = bloom_filter.array_size * CHAR_BIT;
  bloom_filter.hash_num =
    *(table->get_as<decltype(bloom_filter.hash_num)>("hash_num"));
  k = *(table->get_as<decltype(k)>("k"));

  bloom_filter.array = new std::atomic<uint8_t>[bloom_filter.array_size];
  file.read((char*)bloom_filter.array,
            bloom_filter.array_size * sizeof(bloom_filter.array[0]));
}

inline void
KmerBloomFilter::write(const std::string& path)
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
  header->insert("k", get_k());
  root->insert(KMER_BLOOM_FILTER_MAGIC_HEADER, header);
  file << *root << "[HeaderEnd]\n";

  file.write((char*)bloom_filter.array,
             bloom_filter.array_size * sizeof(bloom_filter.array[0]));
}

inline SeedBloomFilter::SeedBloomFilter(size_t bytes,
                                        unsigned k,
                                        const std::vector<std::string>& seeds,
                                        unsigned hash_num_per_seed)
  : kmer_bloom_filter(bytes, hash_num_per_seed, k)
  , seeds(seeds)
  , parsed_seeds(parse_seeds(seeds))
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
    seq, seq_len, get_k(), parsed_seeds, get_hash_num_per_seed());
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
    seq, seq_len, get_k(), parsed_seeds, get_hash_num_per_seed());
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

inline double
SeedBloomFilter::get_fpr() const
{
  const double single_seed_fpr =
    std::pow(get_occupancy(), get_hash_num_per_seed());
  return 1 - std::pow(1 - single_seed_fpr, seeds.size());
}

inline SeedBloomFilter::SeedBloomFilter(const std::string& path)
{
  std::ifstream file(path);

  auto table = kmer_bloom_filter.bloom_filter.parse_header(
    file, SEED_BLOOM_FILTER_MAGIC_HEADER);
  kmer_bloom_filter.bloom_filter.bytes =
    *(table->get_as<decltype(kmer_bloom_filter.bloom_filter.bytes)>("bytes"));
  check_warning(
    sizeof(uint8_t) != sizeof(std::atomic<uint8_t>),
    "Atomic primitives take extra memory. BloomFilter will have less than " +
      std::to_string(get_bytes()) + " for bit array.");
  kmer_bloom_filter.bloom_filter.array_size =
    get_bytes() / sizeof(kmer_bloom_filter.bloom_filter.array[0]);
  kmer_bloom_filter.bloom_filter.array_bits =
    kmer_bloom_filter.bloom_filter.array_size * CHAR_BIT;
  kmer_bloom_filter.bloom_filter.hash_num =
    *(table->get_as<decltype(kmer_bloom_filter.bloom_filter.hash_num)>(
      "hash_num_per_seed"));
  const auto hash_num =
    *(table->get_as<decltype(kmer_bloom_filter.bloom_filter.hash_num)>(
      "hash_num"));
  kmer_bloom_filter.k = *(table->get_as<decltype(kmer_bloom_filter.k)>("k"));
  seeds = *(table->get_array_of<std::string>("seeds"));
  parsed_seeds = parse_seeds(seeds);
  check_error(hash_num != get_hash_num_per_seed() * seeds.size(),
              "SeedBloomFilter: hash_num, hash_num_per_seed, or number of "
              "seeds is wrong.");

  kmer_bloom_filter.bloom_filter.array =
    new std::atomic<uint8_t>[kmer_bloom_filter.bloom_filter.array_size];
  file.read((char*)kmer_bloom_filter.bloom_filter.array,
            kmer_bloom_filter.bloom_filter.array_size *
              sizeof(kmer_bloom_filter.bloom_filter.array[0]));
}

inline void
SeedBloomFilter::write(const std::string& path)
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
  header->insert("hash_num_per_seed", get_hash_num_per_seed());
  header->insert("k", get_k());
  auto seeds_array = cpptoml::make_array();
  for (const auto& seed : seeds) {
    seeds_array->push_back(seed);
  }
  header->insert("seeds", seeds_array);
  root->insert(SEED_BLOOM_FILTER_MAGIC_HEADER, header);
  file << *root << "[HeaderEnd]\n";

  file.write((char*)kmer_bloom_filter.bloom_filter.array,
             kmer_bloom_filter.bloom_filter.array_size *
               sizeof(kmer_bloom_filter.bloom_filter.array[0]));
}

} // namespace btllib

#endif
