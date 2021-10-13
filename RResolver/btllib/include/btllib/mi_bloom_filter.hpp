#ifndef BTLLIB_MI_BLOOM_FILTER_HPP
#define BTLLIB_MI_BLOOM_FILTER_HPP

#include "nthash.hpp"
#include "status.hpp"

#include <sdsl/bit_vector_il.hpp>
#include <sdsl/rank_support.hpp>

#include <algorithm> // std::random_shuffle
#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>
#include <sys/stat.h> // NOLINT
#include <vector>

namespace btllib {

template<typename T>
class MIBloomFilter
{
public:
  static const T MASK = 1 << (sizeof(T) * 8 - 1);
  static const T ANTI_MASK = (T)~MASK;

  static const T STRAND = 1 << (sizeof(T) * 8 - 2);
  static const T ANTI_STRAND = (T)~STRAND;

  static const T ID_MASK = ANTI_STRAND & ANTI_MASK;

  static const unsigned BLOCKSIZE = 512;

  // Calculates the per frame probability of a random match for single value
  static inline double calc_prob_single_frame(double occupancy,
                                              unsigned hash_num,
                                              double freq,
                                              unsigned allowed_misses)
  {
    double prob_total = 0.0;
    for (unsigned i = hash_num - allowed_misses; i <= hash_num; i++) {
      double prob = n_choose_k(hash_num, i);
      prob *= pow(occupancy, i);
      prob *= pow(1.0 - occupancy, hash_num - i);
      prob *= (1.0 - pow(1.0 - freq, i));
      prob_total += prob;
    }
    return prob_total;
  }

  static inline double calc_prob_single(double occupancy, double freq)
  {
    return occupancy * freq;
  }

  /*
   * Returns an a filter size large enough to maintain an occupancy specified
   */
  static size_t calc_optimal_size(size_t entries,
                                  unsigned hash_num,
                                  double occupancy)
  {
    auto non_64_approx_val =
      size_t(-double(entries) * double(hash_num) / log(1.0 - occupancy));
    const int magic = 64;
    return non_64_approx_val + (magic - non_64_approx_val % magic);
  }

  /*
   * Inserts a set of hash values into an sdsl bitvector and returns the number
   * of collisions Thread safe on the bv, though return values will not be the
   * same run to run
   */
  static unsigned insert(sdsl::bit_vector& bv,
                         const uint64_t* hash_values,
                         unsigned hash_num)
  {
    unsigned colli_count = 0;
    for (unsigned i = 0; i < hash_num; ++i) {
      const int magic = 0x3f;
      uint64_t pos = hash_values[i] % bv.size();
      uint64_t* data_index = bv.data() + (pos >> 6); // NOLINT
      uint64_t bit_mask_value = (uint64_t)1 << (pos & magic);
      colli_count +=
        __sync_fetch_and_or(data_index, bit_mask_value) >> (pos & magic) & 1;
    }
    return colli_count;
  }

  // TODO: include allowed miss in header
#pragma pack(1) // to maintain consistent values across platforms
  struct FileHeader
  {
    char magic[8]; // NOLINT
    uint32_t hlen; // header length (including spaced seeds)
    uint64_t size;
    uint32_t nhash;
    uint32_t kmer;
    uint32_t version;
    //		uint8_t allowed_miss;
  };

  /*
   * Constructor using a prebuilt bitvector
   */
  MIBloomFilter<T>(
    unsigned hash_num,
    unsigned kmer_size,
    sdsl::bit_vector& bv,
    const std::vector<std::string>& seeds = std::vector<std::string>(0))
    : m_d_size(0)
    , m_hash_num(hash_num)
    , m_kmer_size(kmer_size)
    , m_sseeds(seeds)
    , m_prob_saturated(0)
  {
    m_bv = sdsl::bit_vector_il<BLOCKSIZE>(bv);
    bv = sdsl::bit_vector();
    if (!seeds.empty()) {
      m_ss_val = parse_seeds(m_sseeds);
      assert(m_sseeds[0].size() == kmer_size);
      for (auto itr = m_sseeds.begin(); itr != m_sseeds.end(); ++itr) {
        // check if spaced seeds are all the same length
        assert(m_kmer_size == itr->size());
      }
    }
    m_rank_support = sdsl::rank_support_il<1>(&m_bv);
    m_d_size = get_pop();
    m_data = new T[m_d_size]();
  }

  MIBloomFilter<T>(const std::string& filter_file_path)
  {
#pragma omp parallel for default(none) shared(filter_file_path)
    for (unsigned i = 0; i < 2; ++i) {
      if (i == 0) {
        FILE* file = fopen(filter_file_path.c_str(), "rbe");
        check_error(file == nullptr,
                    "MIBloomFilter: File " + filter_file_path +
                      " could not be read.");

        FileHeader header;
        check_error(fread(&header, sizeof(struct FileHeader), 1, file) != 1,
                    "MIBloomFilter: Failed to load header.");
        log_info("MIBloomFilter: Loading header...");

        const int magic_nine = 9;
        char magic[magic_nine];
        const int magic_eight = 8;
        memcpy(magic, header.magic, magic_eight);
        magic[magic_eight] = '\0';

        log_info("MIBloomFilter: Loaded header\nmagic: " + std::string(magic) +
                 "\nhlen: " + std::to_string(header.hlen) +
                 "\nsize: " + std::to_string(header.size) +
                 "\nnhash: " + std::to_string(header.nhash) +
                 "\nkmer: " + std::to_string(header.kmer));

        m_hash_num = header.nhash;
        m_kmer_size = header.kmer;
        m_d_size = header.size;
        m_data = new T[m_d_size]();

        if (header.hlen > sizeof(struct FileHeader)) {
          // load seeds
          for (unsigned i = 0; i < header.nhash; ++i) {
            char temp[header.kmer];

            check_error(fread(temp, header.kmer, 1, file) != 1,
                        "MIBloomFilter: Failed to load spaced seed string.");
            log_info("MIBloomFilter: Spaced seed " + std::to_string(i) + ": " +
                     std::string(temp, header.kmer));
            m_sseeds.push_back(std::string(temp, header.kmer));
          }

          m_ss_val = parse_seeds(m_sseeds);
          assert(m_sseeds[0].size() == m_kmer_size);
          for (auto itr = m_sseeds.begin(); itr != m_sseeds.end(); ++itr) {
            // check if spaced seeds are all the same length
            assert(m_kmer_size == itr->size());
          }
        }

        check_error(
          header.hlen != (sizeof(FileHeader) + m_kmer_size * m_sseeds.size()),
          "MIBloomFilter: header length: " + std::to_string(header.hlen) +
            " does not match expected length: " +
            std::to_string(sizeof(FileHeader) + m_kmer_size * m_sseeds.size()) +
            " (likely version mismatch).");

        check_error(strcmp(magic, "MIBLOOMF") != 0,
                    "MIBloomFilter: Bloom filter type does not matc.");

        check_error(header.version != MI_BLOOM_FILTER_VERSION,
                    "MIBloomFilter: Bloom filter version does not match: " +
                      std::to_string(header.version) + " expected " +
                      std::to_string(MI_BLOOM_FILTER_VERSION) + ".");

        log_info("MIBloomFilter: Loading data vector");

        long int l_cur_pos = ftell(file);
        fseek(file, 0, 2);
        size_t file_size = ftell(file) - header.hlen;
        fseek(file, l_cur_pos, 0);

        check_error(file_size != m_d_size * sizeof(T),
                    "MIBloomFilter: " + filter_file_path +
                      " does not match size given by its header. Size: " +
                      std::to_string(file_size) + " vs " +
                      std::to_string(m_d_size * sizeof(T)) + " bytes.");

        size_t count_read = fread(m_data, file_size, 1, file);

        check_error(count_read != 1 && fclose(file) != 0,
                    "MIBloomFilter: File " + filter_file_path +
                      " could not be read.");
      }

      else {
        std::string bv_filename = filter_file_path + ".sdsl";
        log_info("MIBloomFilter: Loading sdsl interleaved bit vector from: " +
                 bv_filename);
        load_from_file(m_bv, bv_filename);
        m_rank_support = sdsl::rank_support_il<1>(&m_bv);
      }
    }

    log_info("MIBloomFilter: Bit vector size: " + std::to_string(m_bv.size()) +
             "\nPopcount: " + std::to_string(get_pop()));
    // TODO: make more streamlined
    m_prob_saturated =
      pow(double(get_pop_saturated()) / double(get_pop()), m_hash_num);
  }

  /*
   * Stores the filter as a binary file to the path specified
   * Stores uncompressed because the random data tends to
   * compress poorly anyway
   */
  void store(std::string const& filter_file_path) const
  {

#pragma omp parallel for default(none) shared(filter_file_path)
    for (unsigned i = 0; i < 2; ++i) {
      if (i == 0) {
        std::ofstream my_file(filter_file_path.c_str(),
                              std::ios::out | std::ios::binary);

        assert(my_file);
        write_header(my_file);

        //				std::cerr << "Storing filter. Filter is
        //"
        //<<
        // m_d_size * sizeof(T)
        //						<< " bytes." <<
        // std::endl;

        // write out each block
        my_file.write(reinterpret_cast<char*>(m_data), m_d_size * sizeof(T));

        my_file.close();
        assert(my_file);

        FILE* file = fopen(filter_file_path.c_str(), "rbe");
        check_error(file == nullptr,
                    "MIBloomFilter: " + filter_file_path +
                      " could not be read.");
      } else {
        std::string bv_filename = filter_file_path + ".sdsl";
        //				std::cerr << "Storing sdsl interleaved
        // bit
        // vector to: " << bv_filename
        //						<< std::endl;
        store_to_file(m_bv, bv_filename);
        //				std::cerr << "Number of bit vector
        // buckets is
        //"
        //<< m_bv.size()
        //						<< std::endl;
        //				std::cerr << "Uncompressed bit vector
        // size is
        //"
        //						<< (m_bv.size() +
        // m_bv.size()
        //* 64
        /// BLOCKSIZE) / 8
        //						<< " bytes" <<
        // std::endl;
      }
    }
  }

  /*
   * Returns false if unable to insert hashes values
   * Contains strand information
   * Inserts hash functions in random order
   */
  bool insert(const uint64_t* hashes, const bool* strand, T val, unsigned max)
  {
    unsigned count = 0;
    std::vector<unsigned> hash_order;
    bool saturated = true;
    // for random number generator seed
    uint64_t rand_value = val;
    bool strand_dir = true;
    if (max % 2 == 0) {
      strand_dir = false;
    }

    // check values and if value set
    for (unsigned i = 0; i < m_hash_num; ++i) {
      // check if values are already set
      uint64_t pos = m_rank_support(hashes[i] % m_bv.size());
      T value = strand_dir ^ strand[i] ? val | STRAND : val;
      // check for saturation
      T old_val = m_data[pos];

      if (old_val > MASK) {
        old_val = old_val & ANTI_MASK;
      } else {
        saturated = false;
      }

      if (old_val == value) {
        ++count;
      } else {
        hash_order.push_back(i);
      }

      if (count >= max) {
        return true;
      }
      rand_value ^= hashes[i];
    }
    std::minstd_rand g(rand_value);
    std::shuffle(hash_order.begin(), hash_order.end(), g);

    // insert seeds in random order
    for (const auto& o : hash_order) {
      uint64_t pos = m_rank_support(hashes[o] % m_bv.size());
      T value = strand_dir ^ strand[o] ? val | STRAND : val;
      // check for saturation
      T old_val = set_val(&m_data[pos], value);

      if (old_val > MASK) {
        old_val = old_val & ANTI_MASK;
      } else {
        saturated = false;
      }

      if (old_val == 0) {
        ++count;
      }

      if (count >= max) {
        return true;
      }
    }

    if (count == 0) {
      if (!saturated) {
        assert(
          max ==
          1); // if this triggers then spaced seed is probably not symmetric
        saturate(hashes);
      }
      return false;
    }
    return true;
  }

  /*
   * Returns false if unable to insert hashes values
   * Inserts hash functions in random order
   */
  bool insert(const uint64_t* hashes, T value, unsigned max)
  {
    unsigned count = 0;
    std::vector<unsigned> hash_order;
    // for random number generator seed
    uint64_t rand_value = value;

    bool saturated = true;

    // check values and if value set
    for (unsigned i = 0; i < m_hash_num; ++i) {
      // check if values are already set
      uint64_t pos = m_rank_support(hashes[i] % m_bv.size());
      // check for saturation
      T old_val = m_data[pos];

      if (old_val > MASK) {
        old_val = old_val & ANTI_MASK;
      } else {
        saturated = false;
      }

      if (old_val == value) {
        ++count;
      } else {
        hash_order.push_back(i);
      }

      if (count >= max) {
        return true;
      }

      rand_value ^= hashes[i];
    }
    std::minstd_rand g(rand_value);
    std::shuffle(hash_order.begin(), hash_order.end(), g);

    // insert seeds in random order
    for (const auto& o : hash_order) {
      uint64_t pos = m_rank_support(hashes[o] % m_bv.size());
      // check for saturation
      T old_val = set_val(&m_data[pos], value);

      if (old_val > MASK) {
        old_val = old_val & ANTI_MASK;
      } else {
        saturated = false;
      }

      if (old_val == 0) {
        ++count;
      }

      if (count >= max) {
        return true;
      }
    }

    if (count == 0) {
      if (!saturated) {
        assert(
          max ==
          1); // if this triggers then spaced seed is probably not symmetric
        saturate(hashes);
      }
      return false;
    }
    return true;
  }

  void saturate(const uint64_t* hashes)
  {
    for (unsigned i = 0; i < m_hash_num; ++i) {
      uint64_t pos = m_rank_support(hashes[i] % m_bv.size());
      __sync_or_and_fetch(&m_data[pos], MASK);
    }
  }

  inline std::vector<T> at(const uint64_t* hashes,
                           bool& saturated,
                           unsigned max_miss = 0)
  {
    std::vector<T> results(m_hash_num);
    unsigned misses = 0;
    for (unsigned i = 0; i < m_hash_num; ++i) {
      uint64_t pos = hashes[i] % m_bv.size();
      if (m_bv[pos] == 0) {
        ++misses;
        saturated = false;
        if (misses > max_miss) {
          return std::vector<T>();
        }
      } else {
        uint64_t rank_pos = m_rank_support(pos);
        T temp_result = m_data[rank_pos];
        if (temp_result > MASK) {
          results[i] = m_data[rank_pos] & ANTI_MASK;
        } else {
          results[i] = m_data[rank_pos];
          saturated = false;
        }
      }
    }
    return results;
  }

  /*
   * Populates rank pos vector. Boolean vector is use to confirm if hits are
   * good Returns total number of misses found
   */
  unsigned at_rank(const uint64_t* hashes,
                   std::vector<uint64_t>& rank_pos,
                   std::vector<bool>& hits,
                   unsigned max_miss) const
  {
    unsigned misses = 0;
    for (unsigned i = 0; i < m_hash_num; ++i) {
      uint64_t pos = hashes[i] % m_bv.size();
      if (bool(m_bv[pos])) {
        rank_pos[i] = m_rank_support(pos);
        hits[i] = true;
      } else {
        if (++misses > max_miss) {
          return misses;
        }

        hits[i] = false;
      }
    }
    return misses;
  }

  /*
   * For k-mers
   * Returns if match succeeded
   */
  bool at_rank(const uint64_t* hashes, std::vector<uint64_t>& rank_pos) const
  {
    for (unsigned i = 0; i < m_hash_num; ++i) {
      uint64_t pos = hashes[i] % m_bv.size();
      if (bool(m_bv[pos])) {
        rank_pos[i] = m_rank_support(pos);
      } else {
        return false;
      }
    }
    return true;
  }

  std::vector<uint64_t> get_rank_pos(const uint64_t* hashes) const
  {
    std::vector<uint64_t> rank_pos(m_hash_num);
    for (unsigned i = 0; i < m_hash_num; ++i) {
      uint64_t pos = hashes[i] % m_bv.size();
      rank_pos[i] = m_rank_support(pos);
    }
    return rank_pos;
  }

  uint64_t get_rank_pos(const uint64_t hash) const
  {
    return m_rank_support(hash % m_bv.size());
  }

  const std::vector<std::vector<unsigned>>& get_seed_values() const
  {
    return m_ss_val;
  }

  unsigned get_kmer_size() const { return m_kmer_size; }

  unsigned get_hash_num() const { return m_hash_num; }

  /*
   * Computes id frequency based on data vector contents
   * Returns counts of repetitive sequence
   */
  size_t get_id_counts(std::vector<size_t>& counts) const
  {
    size_t saturated_counts = 0;
    for (size_t i = 0; i < m_d_size; ++i) {
      if (m_data[i] > MASK) {
        ++counts[m_data[i] & ANTI_MASK];
        ++saturated_counts;
      } else {
        ++counts[m_data[i]];
      }
    }
    return saturated_counts;
  }

  /*
   * computes id frequency based on datavector
   * Returns counts of repetitive sequence
   */
  size_t get_id_counts_strand(std::vector<size_t>& counts) const
  {
    size_t saturated_counts = 0;
    for (size_t i = 0; i < m_d_size; ++i) {
      if (m_data[i] > MASK) {
        ++counts[m_data[i] & ID_MASK];
        ++saturated_counts;
      } else {
        ++counts[m_data[i] & ANTI_STRAND];
      }
    }
    return saturated_counts;
  }

  size_t get_pop() const
  {
    size_t index = m_bv.size() - 1;
    while (m_bv[index] == 0) {
      --index;
    }
    return m_rank_support(index) + 1;
  }

  /*
   * Mostly for debugging
   * should equal get_pop if fully populated
   */
  size_t get_pop_non_zero() const
  {
    size_t count = 0;
    for (size_t i = 0; i < m_d_size; ++i) {
      if (m_data[i] != 0) {
        ++count;
      }
    }
    return count;
  }

  /*
   * Checks data array for abnormal IDs
   * (i.e. values greater than what is possible)
   * Returns first abnormal ID or value of max_val if no abnormal IDs are found
   * For debugging
   */
  T check_values(T max_val) const
  {
    for (size_t i = 0; i < m_d_size; ++i) {
      if ((m_data[i] & ANTI_MASK) > max_val) {
        return m_data[i];
      }
    }
    return max_val;
  }

  size_t get_pop_saturated() const
  {
    size_t count = 0;
    for (size_t i = 0; i < m_d_size; ++i) {
      if (m_data[i] > MASK) {
        ++count;
      }
    }
    return count;
  }

  size_t size() const { return m_bv.size(); }

  // overwrites existing value CAS
  void set_data(uint64_t pos, T id)
  {
    T old_value;
    do {
      old_value = m_data[pos];
      if (old_value > MASK) {
        id |= MASK;
      }
    } while (!__sync_bool_compare_and_swap(&m_data[pos], old_value, id));
  }

  // saturates values
  void saturate_data(uint64_t pos)
  {
#pragma omp critical
    m_data[pos] |= MASK;
  }

  // Does not overwrite
  void set_data_if_empty(uint64_t pos, T id) { set_val(&m_data[pos], id); }

  std::vector<T> get_data(const std::vector<uint64_t>& rank_pos) const
  {
    std::vector<T> results(rank_pos.size());
    for (unsigned i = 0; i < m_hash_num; ++i) {
      results[i] = m_data[rank_pos[i]];
    }
    return results;
  }

  T get_data(uint64_t rank) const { return m_data[rank]; }

  /*
   * Preconditions:
   * 	frame_probs but be equal in size to multiMatchProbs
   * 	frame_probs must be preallocated to correct size (number of ids + 1)
   * Max value is the largest value seen in your set of possible values
   * Returns proportion of saturated elements relative to all elements
   */
  double calc_frame_probs(std::vector<double>& frame_probs,
                          unsigned allowed_miss)
  {
    double occupancy = double(get_pop()) / double(size());
    std::vector<size_t> count_table =
      std::vector<size_t>(frame_probs.size(), 0);
    double sat_prop = double(get_id_counts(count_table));
    size_t sum = 0;
    for (size_t i = 1; i < count_table.size(); ++i) {
      sum += count_table[i];
    }
    sat_prop /= double(sum);
    for (size_t i = 1; i < count_table.size(); ++i) {
      frame_probs[i] =
        calc_prob_single_frame(occupancy,
                               m_hash_num,
                               double(count_table[i]) / double(sum),
                               allowed_miss);
    }
    return sat_prop;
  }

  /*
   * Preconditions:
   * 	frame_probs but be equal in size to multiMatchProbs
   * 	frame_probs must be preallocated to correct size (number of ids + 1)
   * Max value is the largest value seen in your set of possible values
   * Returns proportion of saturated elements relative to all elements
   */
  double calc_frame_probs_strand(std::vector<double>& frame_probs,
                                 unsigned allowed_miss)
  {
    double occupancy = double(get_pop()) / double(size());
    std::vector<size_t> count_table =
      std::vector<size_t>(frame_probs.size(), 0);
    double sat_prop = double(get_id_counts_strand(count_table));
    size_t sum = 0;
    for (const auto& c : count_table) {
      sum += c;
    }
    sat_prop /= double(sum);
#pragma omp parallel for default(none) shared(count_table)
    for (size_t i = 1; i < count_table.size(); ++i) {
      frame_probs[i] =
        calc_prob_single_frame(occupancy,
                               m_hash_num,
                               double(count_table[i]) / double(sum),
                               allowed_miss);
      //			frame_probs[i] = calc_prob_single(occupancy,
      //					double(count_table[i]) /
      // double(sum));
    }
    return sat_prop;
  }

  ~MIBloomFilter() { delete[] m_data; }

private:
  // Driver function to sort the std::vector elements
  // by second element of pairs
  static bool sort_by_sec(const std::pair<int, int>& a,
                          const std::pair<int, int>& b)
  {
    return (a.second < b.second);
  }

  /*
   * Helper function for header storage
   */
  void write_header(std::ofstream& out) const
  {
    FileHeader header;
    const int magic_num = 8;
    memcpy(header.magic, "MIBLOOMF", magic_num);

    header.hlen = sizeof(struct FileHeader) + m_kmer_size * m_sseeds.size();
    header.kmer = m_kmer_size;
    header.size = m_d_size;
    header.nhash = m_hash_num;
    header.version = MI_BLOOM_FILTER_VERSION;

    //		std::cerr << "Writing header... magic: " << magic << " hlen: "
    //<<
    // header.hlen
    //				<< " nhash: " << header.nhash << " size: " <<
    // header.size
    //				<< std::endl;

    out.write(reinterpret_cast<char*>(&header), sizeof(struct FileHeader));

    for (const auto& s : m_sseeds) {
      out.write(s.c_str(), m_kmer_size);
    }
  }

  /*
   * Calculates the optimal number of hash function to use
   * Calculation assumes optimal ratio of bytes per entry given a fpr
   */
  inline static unsigned calc_opti_hash_num(double fpr)
  {
    return unsigned(-log(fpr) / log(2));
  }

  /*
   * Calculate FPR based on hash functions, size and number of entries
   * see http://en.wikipedia.org/wiki/Bloom_filter
   */
  double calc_fpr_num_inserted(size_t num_entr) const
  {
    return pow(1.0 - pow(1.0 - 1.0 / double(m_bv.size()),
                         double(num_entr) * double(m_hash_num)),
               double(m_hash_num));
  }

  /*
   * Calculates the optimal FPR to use based on hash functions
   */
  double calc_fpr_hash_num(int hash_funct_num) const
  {
    const double magic = 2.0;
    return pow(magic, -hash_funct_num);
  }

  /*
   * Returns old value that was inside
   * Does not overwrite if non-zero value already exists
   */
  T set_val(T* val, T new_val)
  {
    T old_value;
    do {
      old_value = *val;
      if (old_value != 0) {
        break;
      }
    } while (!__sync_bool_compare_and_swap(val, old_value, new_val));
    return old_value;
  }

  static inline unsigned n_choose_k(unsigned n, unsigned k)
  {
    if (k > n) {
      return 0;
    }
    if (k * 2 > n) {
      k = n - k;
    }
    if (k == 0) {
      return 1;
    }
    unsigned result = n;
    for (unsigned i = 2; i <= k; ++i) {
      result *= (n - i + 1);
      result /= i;
    }
    return result;
  }

  // size of bitvector
  size_t m_d_size;

  sdsl::bit_vector_il<BLOCKSIZE> m_bv;
  T* m_data;
  sdsl::rank_support_il<1> m_rank_support;

  unsigned m_hash_num;
  unsigned m_kmer_size;

  using seed_val = std::vector<std::vector<unsigned>>;
  std::vector<std::string> m_sseeds;

  double m_prob_saturated;
  seed_val m_ss_val;

  static const uint32_t MI_BLOOM_FILTER_VERSION = 1;
};

} // namespace btllib

#endif
