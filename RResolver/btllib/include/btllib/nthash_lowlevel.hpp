/*
 * nthash_lowlevel.hpp
 * Author: Hamid Mohamadi
 * Genome Sciences Centre,
 * British Columbia Cancer Agency
 */

#ifndef BTLLIB_NTHASH_LOWLEVEL_HPP
#define BTLLIB_NTHASH_LOWLEVEL_HPP

#include "btllib/nthash_consts.hpp"
#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

namespace btllib {

template<typename T>
inline T
canonical(const T fwd, const T rev)
{
  return fwd + rev;
}

static_assert(std::numeric_limits<unsigned>::max() + 1 == 0,
              "Integers don't overflow on this platform which is necessary for "
              "ntHash canonical hash computation.");

// Data structures for spaced seeds

// List of don't care positions
using SpacedSeed = std::vector<unsigned>;

// Bx2 array representing block start and end positions.
using SpacedSeedBlocks = std::vector<std::array<unsigned, 2>>;

// List of blocks with sizes of one.
using SpacedSeedMonomers = std::vector<unsigned>;

/**
 * Split a 64-bit word into 33 and 31-bit subwords and left-rotate them
 * separately.
 *
 * @param x A 64-bit unsigned integer.
 *
 * @return Split-rotation result.
 */
inline uint64_t
srol(const uint64_t x)
{
  uint64_t m = ((x & 0x8000000000000000ULL) >> 30) | // NOLINT
               ((x & 0x100000000ULL) >> 32);         // NOLINT
  return ((x << 1) & 0xFFFFFFFDFFFFFFFFULL) | m;     // NOLINT
}

/**
 * Split a 64-bit word into 33 and 31-bit subwords and left-rotate them
 * separately multiple times.
 *
 * @param x A 64-bit unsigned integer.
 * @param d Number of rotations.
 *
 * @return Split-rotation result.
 */
inline uint64_t
srol(const uint64_t x, const unsigned d)
{
  uint64_t v = (x << d) | (x >> (64 - d));                         // NOLINT
  uint64_t y = (v ^ (v >> 33)) &                                   // NOLINT
               (std::numeric_limits<uint64_t>::max() >> (64 - d)); // NOLINT
  return v ^ (y | (y << 33));                                      // NOLINT
}

/**
 * Split a 64-bit word into 33 and 31-bit subwords and right-rotate them
 * separately.
 *
 * @param x A 64-bit unsigned integer.
 *
 * @return Split-rotation result.
 */
inline uint64_t
sror(const uint64_t x)
{
  uint64_t m = ((x & 0x200000000ULL) << 30) | ((x & 1ULL) << 32); // NOLINT
  return ((x >> 1) & 0xFFFFFFFEFFFFFFFFULL) | m;                  // NOLINT
}

/**
 * Generate the forward-strand hash value of the first k-mer in the sequence.
 *
 * @param kmer_seq C array containing the sequence's characters.
 * @param k k-mer size.
 *
 * @return Hash value of k-mer_0.
 */
uint64_t
ntf64(const char* kmer_seq, unsigned k);

/**
 * Generate a hash value for the reverse-complement of the first k-mer in the
 * sequence.
 *
 * @param kmer_seq C array containing the sequence's characters.
 * @param k k-mer size.
 *
 * @return Hash value of the reverse-complement of k-mer_0.
 */
uint64_t
ntr64(const char* kmer_seq, unsigned k);

/**
 * Perform a roll operation on the forward strand by removing char_out and
 * including char_in.
 *
 * @param fh_val Previous hash value computed for the sequence.
 * @param k k-mer size.
 * @param char_out Character to be removed.
 * @param char_in Character to be included.
 *
 * @return Rolled forward hash value.
 */
uint64_t
ntf64(uint64_t fh_val,
      unsigned k,
      unsigned char char_out,
      unsigned char char_in);

/**
 * Perform a roll operation on the reverse-complement by removing char_out and
 * including char_in.
 *
 * @param rh_val Previous reverse-complement hash value computed for the
 * sequence.
 * @param k k-mer size.
 * @param char_out Character to be removed.
 * @param char_in Character to be included.
 *
 * @return Rolled hash value for the reverse-complement.
 */
uint64_t
ntr64(uint64_t rh_val,
      unsigned k,
      unsigned char char_out,
      unsigned char char_in);

/**
 * Generate a canonical hash value for the first k-mer.
 *
 * @param kmer_seq C array containing the sequence's characters.
 * @param k k-mer size.
 *
 * @return Canonical hash value of k-mer_0.
 */
uint64_t
ntc64(const char* kmer_seq, unsigned k);

/**
 * Generate a canonical hash value for the first k-mer and update both strands'
 * hash values.
 *
 * @param kmer_seq C array containing the sequence's characters.
 * @param k k-mer size.
 * @param fh_val Forward strand hash value container.
 * @param rh_val Reverse strand hash value container.
 *
 * @return Canonical hash value of k-mer_0.
 */
uint64_t
ntc64(const char* kmer_seq, unsigned k, uint64_t& fh_val, uint64_t& rh_val);

/**
 * Perform a roll operation on the sequence and generate a canonical hash value.
 *
 * @param char_out Character to be removed.
 * @param char_in Character to be included.
 * @param k k-mer size.
 * @param fh_val Previous hash value for the forward strand.
 * @param rh_val Previous hash value for the reverse-complement.
 *
 * @return Canonical hash value after including char_in and removing char_out.
 */
uint64_t
ntc64(unsigned char char_out,
      unsigned char char_in,
      unsigned k,
      uint64_t& fh_val,
      uint64_t& rh_val);

/**
 * Perform a roll-back operation on the forward strand.
 *
 * @param rh_val Previous forward hash value computed for the sequence.
 * @param k k-mer size.
 * @param char_out Character to be removed.
 * @param char_in Character to be included.
 *
 * @return Resulting hash value.
 */
uint64_t
ntf64l(uint64_t rh_val,
       unsigned k,
       unsigned char char_out,
       unsigned char char_in);

/**
 * Perform a roll-back operation on the reverse-complement.
 *
 * @param rh_val Previous reverse hash value computed for the sequence.
 * @param k k-mer size.
 * @param char_out Character to be removed.
 * @param char_in Character to be included.
 *
 * @return Resulting hash value for the reverse-complement.
 */
uint64_t
ntr64l(uint64_t fh_val,
       unsigned k,
       unsigned char char_out,
       unsigned char char_in);

/**
 * Perform a roll-back operation on the canonical hash value and update previous
 * hashes for both strands.
 *
 * @param char_out Character to be removed.
 * @param char_in Character to be included.
 * @param k k-mer size.
 * @param fh_val Previous forward hash value computed for the sequence.
 * @param rh_val Previous reverse hash value computed for the sequence.
 *
 * @return Roll back result for the canonical hash value.
 */
uint64_t
ntc64l(unsigned char char_out,
       unsigned char char_in,
       unsigned k,
       uint64_t& fh_val,
       uint64_t& rh_val);

/**
 * Extend hash array using a base hash value.
 *
 * @param bh_val Base hash value.
 * @param k k-mer size.
 * @param h Size of the resulting hash array (number of extra hashes minus one).
 * @param h_val Array of size h for storing the output hashes.
 */
void
nte64(uint64_t bh_val, unsigned k, unsigned h, uint64_t* h_val);

/**
 * Generate multiple canonical hash values for the first k-mer.
 *
 * @param kmer_seq Array containing the sequence's characters.
 * @param k k-mer size.
 * @param m Number of hashes per k-mer.
 * @param h_val Array of size m for storing the hash values.
 */
void
ntmc64(const char* kmer_seq, unsigned k, unsigned m, uint64_t* h_val);

/**
 * Generate multiple canonical hash values for the first k-mer and return
 * strand-specific hash values.
 *
 * @param kmer_seq Array containing the sequence's characters.
 * @param k k-mer size.
 * @param m Number of hashes per k-mer.
 * @param fh_val Unsigned 64-bit int container for the forward hash.
 * @param rh_val Unsigned 64-bit int container for the reverse-complement hash.
 * @param h_val Array of size m for storing the hash values.
 */
void
ntmc64(const char* kmer_seq,
       unsigned k,
       unsigned m,
       uint64_t& fh_val,
       uint64_t& rh_val,
       uint64_t* h_val);

/**
 * Generate a new canonical hash value by performing a roll operation.
 *
 * @param char_out Character to be removed.
 * @param char_in Character to be included.
 * @param k k-mer size.
 * @param m Number of hashes per k-mer.
 * @param fh_val Previous forward hash value.
 * @param rh_val Previous reverse hash value.
 * @param h_val Array of size m for storing the output hash values.
 */
void
ntmc64(unsigned char char_out,
       unsigned char char_in,
       unsigned k,
       unsigned m,
       uint64_t& fh_val,
       uint64_t& rh_val,
       uint64_t* h_val);

/**
 * Generate a new canonical hash value by performing a roll-back operation.
 *
 * @param char_out Character to be removed.
 * @param char_in Character to be included.
 * @param k k-mer size.
 * @param m Number of hashes per k-mer.
 * @param fh_val Previous forward hash value.
 * @param rh_val Previous reverse hash value.
 * @param h_val Array of size m for storing the output hash values.
 */
void
ntmc64l(unsigned char char_out,
        unsigned char char_in,
        unsigned k,
        unsigned m,
        uint64_t& fh_val,
        uint64_t& rh_val,
        uint64_t* h_val);

/**
 * Generate a canonical hash value for the first k-mer and find the first
 * ignored character.
 *
 * @param kmer_seq Array containing the sequence's characters.
 * @param k k-mer size.
 * @param h_val Container for the output hash value.
 * @param loc_n Location of the first unknown character.
 *
 * @return true if all the characters of the first k-mer are known, otherwise
 * false.
 */
bool
ntc64(const char* kmer_seq, unsigned k, uint64_t& h_val, unsigned& loc_n);

/**
 * Generate multiple canonical hash values for the first k-mer and find the
 * first ignored character.
 *
 * @param kmer_seq Array containing the sequence's characters.
 * @param k k-mer size.
 * @param m Number of hashes per k-mer.
 * @param h_val Array of size m for storing the output hash values.
 * @param loc_n Location of the first unknown character.
 *
 * @return true if all the characters of the first k-mer are known, otherwise
 * false.
 */
bool
ntmc64(const char* kmer_seq,
       unsigned k,
       unsigned m,
       unsigned& loc_n,
       uint64_t* h_val);

/**
 * Generate a canonical hash value for the first k-mer, find the first ignored
 * character and return the strand-specific hash values.
 *
 * @param kmer_seq Array containing the sequence's characters.
 * @param k k-mer size.
 * @param fh_val Container for the forward hash value.
 * @param rh_val Container for the reverse hash value.
 * @param h_val Container for the output hash value.
 * @param loc_n Location of the first unknown character.
 *
 * @return true if all the characters of the first k-mer are known, otherwise
 * false.
 */
bool
ntc64(const char* kmer_seq,
      unsigned k,
      uint64_t& fh_val,
      uint64_t& rh_val,
      uint64_t& h_val,
      unsigned& loc_n);

/**
 * Generate multiple canonical hash value for the first k-mer, find the first
 * ignored character and return the strand-specific hash values.
 *
 * @param kmer_seq Array containing the sequence's characters.
 * @param k k-mer size.
 * @param m Number of hashes per k-mer.
 * @param fh_val Container for the forward hash value.
 * @param rh_val Container for the reverse hash value.
 * @param loc_n Location of the first unknown character.
 * @param h_val Array of size m for storing the output hash values.
 *
 * @return true if all the characters of the first k-mer are known, otherwise
 * false.
 */
bool
ntmc64(const char* kmer_seq,
       unsigned k,
       unsigned m,
       uint64_t& fh_val,
       uint64_t& rh_val,
       unsigned& loc_n,
       uint64_t* h_val);

/**
 * Generate multiple canonical hash values for the first k-mer, find the first
 * ignored character, and returning the strand-specific hash values and strand
 * selections.
 *
 * @param kmer_seq Array containing the sequence's characters.
 * @param k k-mer size.
 * @param m Number of hashes per k-mer.
 * @param fh_val Container for the forward hash value.
 * @param rh_val Container for the reverse hash value.
 * @param loc_n Location of the first unknown character.
 * @param h_val Array of size m for storing the output hash values.
 * @param h_stn true if the reverse strand was selected, otherwise false.
 *
 * @return true if all the characters of the first k-mer are known, otherwise
 * false.
 */
bool
ntmc64(const char* kmer_seq,
       unsigned k,
       unsigned m,
       uint64_t& fh_val,
       uint64_t& rh_val,
       unsigned& loc_n,
       uint64_t* h_val,
       bool& h_stn);

/**
 * Generate multiple canonical hash values by performing a roll operation,
 * returning the strand-specific hash values and strand selections.
 *
 * @param char_out Character to be removed.
 * @param char_in Character to be included.
 * @param k k-mer size.
 * @param m Number of hashes per k-mer.
 * @param fh_val Container for the forward hash value.
 * @param rh_val Container for the reverse hash value.
 * @param h_val Array of size m for storing the output hash values.
 * @param h_stn true if the reverse strand was selected, otherwise false.
 */
void
ntmc64(unsigned char char_out,
       unsigned char char_in,
       unsigned k,
       unsigned m,
       uint64_t& fh_val,
       uint64_t& rh_val,
       uint64_t* h_val,
       bool& h_stn);

/**
 * Generate a hash value for the input spaced seed by excluding all don't care
 * positions.
 *
 * @param fk_val Forward hash value of the k-mer (ignoring the spaced seed).
 * @param rk_val Reverse hash value of the k-mer (ignoring the spaced seed).
 * @param seed_seq Array of characters representing the spaced seed. Anything
 * other than '1' is treated as a don't care.
 * @param kmer_seq Array of character representing the k-mer.
 * @param k k-mer size.
 *
 * @return Canonical hash value for the k-mer masked with the spaced seed.
 */
uint64_t
mask_hash(uint64_t& fk_val,
          uint64_t& rk_val,
          const char* seed_seq,
          const char* kmer_seq,
          unsigned k);

/**
 * Generate multiple new hash values for the input k-mer by substituting
 * multiple characters.
 *
 * @param fh_val Forward hash value of the k-mer.
 * @param rh_val Reverse hash value of the k-mer.
 * @param kmer_seq Array of characters representing the k-mer.
 * @param positions Indicies of the positions to be substituted.
 * @param new_bases Characters to be placed in the indicies indicated in
 * positions.
 * @param k k-mer size.
 * @param m Number of hashes per k-mer.
 * @param h_val Array of size m for storing the output hash values.
 */
void
sub_hash(uint64_t fh_val,
         uint64_t rh_val,
         const char* kmer_seq,
         const std::vector<unsigned>& positions,
         const std::vector<unsigned char>& new_bases,
         unsigned k,
         unsigned m,
         uint64_t* h_val);

/**
 * Generate multiple hash values for the input spaced seeds and first k-mer.
 *
 * @param kmer_seq Array of characters representing the k-mer.
 * @param seed_seq Array of SpacedSeed objects representing the seeds' blocks.
 * @param monomers List of the positions that represent blocks of size one for
 * each seed.
 * @param k k-mer size.
 * @param m Number of spaced seeds.
 * @param m2 Number of hashes per seed.
 * @param fh_nomonos Container for the forward hash values before including the
 * size-one blocks.
 * @param rh_nomonos Container for the reverse hash values before including the
 * size-one blocks.
 * @param fh_val Container for the forward hash values after including the
 * size-one blocks.
 * @param rh_val Container for the reverse hash values after including the
 * size-one blocks.
 * @param loc_n Location of the first unknown character in the first sequence.
 * @param h_val Array of size m * m2 for storing the output hash values.
 *
 * @return true if all the care positions of the first k-mer are valid,
 * otherwise false.
 */
bool
ntmsm64(const char* kmer_seq,
        const std::vector<SpacedSeedBlocks>& seeds_blocks,
        const std::vector<SpacedSeedMonomers>& seeds_monomers,
        unsigned k,
        unsigned m,
        unsigned m2,
        uint64_t* fh_nomonos,
        uint64_t* rh_nomonos,
        uint64_t* fh_val,
        uint64_t* rh_val,
        unsigned& loc_n,
        uint64_t* h_val);

#define NTMSM64(ROL_HANDLING, IN_HANDLING, OUT_HANDLING, ROR_HANDLING)         \
  unsigned char char_out, char_in;                                             \
  uint64_t fh_seed, rh_seed;                                                   \
  unsigned i_out, i_in, i_base;                                                \
  for (unsigned i_seed = 0; i_seed < m; i_seed++) {                            \
    ROL_HANDLING /* NOLINT(bugprone-macro-parentheses) */                      \
      for (const auto& block : seeds_blocks[i_seed])                           \
    {                                                                          \
      IN_HANDLING                                                              \
      OUT_HANDLING                                                             \
      fh_seed ^= MS_TAB(char_out, k - i_out);                                  \
      fh_seed ^= MS_TAB(char_in, k - i_in);                                    \
      rh_seed ^= MS_TAB(char_out & CP_OFF, i_out);                             \
      rh_seed ^= MS_TAB(char_in & CP_OFF, i_in);                               \
    }                                                                          \
    ROR_HANDLING /* NOLINT(bugprone-macro-parentheses) */                      \
      fh_nomonos[i_seed] = fh_seed;                                            \
    rh_nomonos[i_seed] = rh_seed;                                              \
    for (const auto& pos : seeds_monomers[i_seed]) {                           \
      fh_seed ^= MS_TAB((unsigned char)kmer_seq[pos + 1], k - 1 - pos);        \
      rh_seed ^= MS_TAB((unsigned char)kmer_seq[pos + 1] & CP_OFF, pos);       \
    }                                                                          \
    fh_val[i_seed] = fh_seed;                                                  \
    rh_val[i_seed] = rh_seed;                                                  \
    i_base = i_seed * m2;                                                      \
    h_val[i_base] = canonical(fh_seed, rh_seed);                               \
    for (unsigned i_hash = 1; i_hash < m2; i_hash++) {                         \
      h_val[i_base + i_hash] = h_val[i_base] * (i_hash ^ k * MULTISEED);       \
      h_val[i_base + i_hash] ^= h_val[i_base + i_hash] >> MULTISHIFT;          \
    }                                                                          \
  }

/**
 * Generate multiple hash values for the input spaced seeds and the next
 * k-mer by performing a forward roll operation.
 *
 * @param kmer_seq Array of characters representing the previous k-mer.
 * @param seed_seq Array of SpacedSeed objects representing the seeds' blocks.
 * @param monomers List of the positions that represent blocks of size one for
 * each seed.
 * @param k k-mer size.
 * @param m Number of spaced seeds.
 * @param m2 Number of hashes per seed.
 * @param fh_nomonos Previous forward hash values before including the size-one
 * blocks.
 * @param rh_nomonos Previous reverse hash values before including the size-one
 * blocks.
 * @param fh_val Previous forward hash values after including the size-one
 * blocks.
 * @param rh_val Previous reverse hash values after including the size-one
 * blocks.
 * @param h_val Array of size m * m2 for storing the output hash values.
 */
void
ntmsm64(const char* kmer_seq,
        const std::vector<SpacedSeedBlocks>& seeds_blocks,
        const std::vector<SpacedSeedMonomers>& seeds_monomers,
        unsigned k,
        unsigned m,
        unsigned m2,
        uint64_t* fh_nomonos,
        uint64_t* rh_nomonos,
        uint64_t* fh_val,
        uint64_t* rh_val,
        uint64_t* h_val);

/**
 * Generate multiple hash values for the input spaced seeds and the next
 * k-mer by performing a backward roll operation.
 *
 * @param kmer_seq Array of characters representing the previous k-mer.
 * @param seed_seq Array of SpacedSeed objects representing the seeds' blocks.
 * @param monomers List of the positions that represent blocks of size one for
 * each seed.
 * @param k k-mer size.
 * @param m Number of spaced seeds.
 * @param m2 Number of hashes per seed.
 * @param fh_nomonos Previous forward hash values before including the size-one
 * blocks.
 * @param rh_nomonos Previous reverse hash values before including the size-one
 * blocks.
 * @param fh_val Previous forward hash values after including the size-one
 * blocks.
 * @param rh_val Previous reverse hash values after including the size-one
 * blocks.
 * @param h_val Array of size m * m2 for storing the output hash values.
 */
void
ntmsm64l(const char* kmer_seq,
         const std::vector<SpacedSeedBlocks>& seeds_blocks,
         const std::vector<SpacedSeedMonomers>& seeds_monomers,
         unsigned k,
         unsigned m,
         unsigned m2,
         uint64_t* fh_nomonos,
         uint64_t* rh_nomonos,
         uint64_t* fh_val,
         uint64_t* rh_val,
         uint64_t* h_val);

/**
 * Generate multiple hash values for the input spaced seeds and the next
 * k-mer by performing a forward peek operation.
 *
 * @param kmer_seq Array of characters representing the previous k-mer.
 * @param seed_seq Array of SpacedSeed objects representing the seeds' blocks.
 * @param monomers List of the positions that represent blocks of size one for
 * each seed.
 * @param k k-mer size.
 * @param m Number of spaced seeds.
 * @param m2 Number of hashes per seed.
 * @param fh_nomonos Previous forward hash values before including the size-one
 * blocks.
 * @param rh_nomonos Previous reverse hash values before including the size-one
 * blocks.
 * @param fh_val Previous forward hash values after including the size-one
 * blocks.
 * @param rh_val Previous reverse hash values after including the size-one
 * blocks.
 * @param h_val Array of size m * m2 for storing the output hash values.
 */
void
ntmsm64(const char* kmer_seq,
        char in,
        const std::vector<SpacedSeedBlocks>& seeds_blocks,
        const std::vector<SpacedSeedMonomers>& seeds_monomers,
        unsigned k,
        unsigned m,
        unsigned m2,
        uint64_t* fh_nomonos,
        uint64_t* rh_nomonos,
        uint64_t* fh_val,
        uint64_t* rh_val,
        uint64_t* h_val);

/**
 * Generate multiple hash values for the input spaced seeds and the next
 * k-mer by performing a backwards peek operation.
 *
 * @param kmer_seq Array of characters representing the previous k-mer.
 * @param seed_seq Array of SpacedSeed objects representing the seeds' blocks.
 * @param monomers List of the positions that represent blocks of size one for
 * each seed.
 * @param k k-mer size.
 * @param m Number of spaced seeds.
 * @param m2 Number of hashes per seed.
 * @param fh_nomonos Previous forward hash values before including the size-one
 * blocks.
 * @param rh_nomonos Previous reverse hash values before including the size-one
 * blocks.
 * @param fh_val Previous forward hash values after including the size-one
 * blocks.
 * @param rh_val Previous reverse hash values after including the size-one
 * blocks.
 * @param h_val Array of size m * m2 for storing the output hash values.
 */
void
ntmsm64l(const char* kmer_seq,
         char in,
         const std::vector<SpacedSeedBlocks>& seeds_blocks,
         const std::vector<SpacedSeedMonomers>& seeds_monomers,
         unsigned k,
         unsigned m,
         unsigned m2,
         uint64_t* fh_nomonos,
         uint64_t* rh_nomonos,
         uint64_t* fh_val,
         uint64_t* rh_val,
         uint64_t* h_val);

} // namespace btllib

#endif