#ifndef BTLLIB_BLOOM_FILTER_HPP
#define BTLLIB_BLOOM_FILTER_HPP

#include "nthash.hpp"
#include "status.hpp"

#include "vendor/cpptoml.hpp"

#include <atomic>
#include <climits>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <memory>
#include <string>

namespace btllib {

static const unsigned char BIT_MASKS[CHAR_BIT] = {
  // NOLINT
  0x01, 0x02, 0x04, 0x08, // NOLINT
  0x10, 0x20, 0x40, 0x80  // NOLINT
};

static const char* const BLOOM_FILTER_MAGIC_HEADER = "BTLBloomFilter_v2";
static const char* const KMER_BLOOM_FILTER_MAGIC_HEADER =
  "BTLKmerBloomFilter_v2";

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
  BloomFilter() {}
  BloomFilter(size_t bytes, unsigned hash_num);
  BloomFilter(const std::string& path);

  virtual ~BloomFilter() { delete[] array; }

  void insert(const std::vector<uint64_t>& hashes);
  void insert(const uint64_t* hashes);

  bool contains(const std::vector<uint64_t>& hashes) const;
  bool contains(const uint64_t* hashes) const;

  size_t get_bytes() const { return bytes; }
  uint64_t get_pop_cnt() const;
  unsigned get_hash_num() const { return hash_num; }
  double get_fpr() const;

  static std::shared_ptr<cpptoml::table> parse_header(
    std::ifstream& file,
    const std::string& magic_string);

  void write(const std::string& path);

protected:
  std::atomic<uint8_t>* array = nullptr;
  size_t bytes = 0;
  size_t array_size =
    0; // Should be equal to bytes, but not guaranteed by standard
  size_t array_bits = 0;
  unsigned hash_num = 0;
};

/**
 * Bloom filter data structure that kmerizes and hashes given sequences,
 * storing the results.
 */
class KmerBloomFilter : public BloomFilter
{

public:
  /**
   * Constructor.
   * @param k kmer size
   * @param bytes bytes to allocate for the filter
   * @param hash_num number of hashes
   */
  KmerBloomFilter(size_t bytes, unsigned hash_num, unsigned k);
  KmerBloomFilter(const std::string& path);

  ~KmerBloomFilter() override {}

  /**
   * Store the kmers of a sequence.
   * @param seq sequence to kmerize
   */
  void insert(const std::string& seq);

  /**
   * Store the kmers of a sequence.
   * @param seq sequence to kmerize
   * @param seq_len length of seq
   */
  void insert(const char* seq, size_t seq_len);

  /**
   * Query the kmers of a sequence.
   * @param seq sequence to kmerize
   *
   * @return number of kmers found in seq
   */
  unsigned contains(const std::string& seq) const;

  /**
   * Query the kmers of a sequence.
   * @param seq sequence to kmerize
   * @param seq_len length of seq
   *
   * @return number of kmers found in seq
   */
  unsigned contains(const char* seq, size_t seq_len) const;

  void write(const std::string& path);

protected:
  unsigned k;
};

inline BloomFilter::BloomFilter(size_t bytes, unsigned hash_num)
  : bytes(std::ceil(bytes / sizeof(uint64_t)) * sizeof(uint64_t))
  , array_size(this->bytes / sizeof(array[0]))
  , array_bits(array_size * CHAR_BIT)
  , hash_num(hash_num)
{
  check_warning(
    sizeof(uint8_t) != sizeof(std::atomic<uint8_t>),
    "Atomic primitives take extra memory. BloomFilter will have less than " +
      std::to_string(bytes) + " for bit array.");
  array = new std::atomic<uint8_t>[array_size];
  std::memset((void*)array, 0, array_size * sizeof(array[0]));
}

inline void
BloomFilter::insert(const std::vector<uint64_t>& hashes)
{
  insert(hashes.data());
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
BloomFilter::contains(const std::vector<uint64_t>& hashes) const
{
  return contains(hashes.data());
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
#pragma omp parallel for reduction(+ : pop_cnt)
  for (size_t i = 0; i < array_size; ++i) {
    pop_cnt += pop_cnt_byte(array[i]);
  }
  return pop_cnt;
}

inline double
BloomFilter::get_fpr() const
{
  return std::pow(double(get_pop_cnt()) / double(array_size), double(hash_num));
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
  bytes = *table->get_as<decltype(bytes)>("bytes");
  check_warning(
    sizeof(uint8_t) != sizeof(std::atomic<uint8_t>),
    "Atomic primitives take extra memory. BloomFilter will have less than " +
      std::to_string(bytes) + " for bit array.");
  array_size = bytes / sizeof(std::atomic<uint8_t>);
  array_bits = array_size * CHAR_BIT;
  hash_num = *table->get_as<decltype(hash_num)>("hash_num");

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
  header->insert("bytes", bytes);
  header->insert("hash_num", hash_num);
  root->insert(BLOOM_FILTER_MAGIC_HEADER, header);
  file << *root << "[HeaderEnd]\n";

  file.write((char*)array, array_size * sizeof(array[0]));
}

inline KmerBloomFilter::KmerBloomFilter(size_t bytes,
                                        unsigned hash_num,
                                        unsigned k)
  : BloomFilter(bytes, hash_num)
  , k(k)
{}

inline void
KmerBloomFilter::insert(const std::string& seq)
{
  insert(seq.c_str(), seq.size());
}

inline void
KmerBloomFilter::insert(const char* seq, size_t seq_len)
{
  NtHash nthash(seq, seq_len, k, get_hash_num());
  while (nthash.roll()) {
    BloomFilter::insert(nthash.hashes());
  }
}

inline unsigned
KmerBloomFilter::contains(const std::string& seq) const
{
  return contains(seq.c_str(), seq.size());
}

inline unsigned
KmerBloomFilter::contains(const char* seq, size_t seq_len) const
{
  unsigned count = 0;
  NtHash nthash(seq, seq_len, k, get_hash_num());
  while (nthash.roll()) {
    if (BloomFilter::contains(nthash.hashes())) {
      count++;
    }
  }
  return count;
}

inline KmerBloomFilter::KmerBloomFilter(const std::string& path)
{
  std::ifstream file(path);

  auto table = parse_header(file, KMER_BLOOM_FILTER_MAGIC_HEADER);
  bytes = *table->get_as<decltype(bytes)>("bytes");
  check_warning(
    sizeof(uint8_t) != sizeof(std::atomic<uint8_t>),
    "Atomic primitives take extra memory. BloomFilter will have less than " +
      std::to_string(bytes) + " for bit array.");
  array_size = bytes / sizeof(array[0]);
  array_bits = array_size * CHAR_BIT;
  hash_num = *table->get_as<decltype(hash_num)>("hash_num");
  k = *table->get_as<decltype(k)>("k");

  array = new std::atomic<uint8_t>[array_size];
  file.read((char*)array, array_size * sizeof(array[0]));
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
  header->insert("bytes", bytes);
  header->insert("hash_num", hash_num);
  header->insert("k", hash_num);
  root->insert(KMER_BLOOM_FILTER_MAGIC_HEADER, header);
  file << *root << "[HeaderEnd]\n";

  file.write((char*)array, array_size * sizeof(array[0]));
}

} // namespace btllib

#endif
