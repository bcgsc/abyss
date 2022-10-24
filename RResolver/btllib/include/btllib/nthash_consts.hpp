/*
 * nthash_consts.hpp
 * Author: Hamid Mohamadi
 * Genome Sciences Centre,
 * British Columbia Cancer Agency
 */

#ifndef BTLLIB_NTHASH_CONSTS_HPP
#define BTLLIB_NTHASH_CONSTS_HPP

#include <cmath>
#include <cstdint>

namespace btllib {

#define MS_TAB(CHAR, ROT)                                                      \
  (MS_TAB_31L[CHAR][(ROT) < 31 ? (ROT) : (ROT) % 31] | /* NOLINT */            \
   MS_TAB_33R[CHAR][(ROT) < 33 ? (ROT) : (ROT) % 33])  /* NOLINT */

// offset for the complement base in the random seeds table
const uint8_t CP_OFF = 0x07;

// shift for gerenerating multiple hash values
const int MULTISHIFT = 27;

// seed for gerenerating multiple hash values
static const uint64_t MULTISEED = 0x90b45d39fb6da1fa;

// 64-bit random seeds corresponding to bases and their complements
static const uint64_t SEED_A = 0x3c8bfbb395c60474;
static const uint64_t SEED_C = 0x3193c18562a02b4c;
static const uint64_t SEED_G = 0x20323ed082572324;
static const uint64_t SEED_T = 0x295549f54be24456;
static const uint64_t SEED_N = 0x0000000000000000;

static const int ASCII_SIZE = 256;

extern const uint64_t SEED_TAB[ASCII_SIZE];

extern const uint64_t A33R[33];
extern const uint64_t A31L[31];

extern const uint64_t C33R[33];
extern const uint64_t C31L[31];

extern const uint64_t G33R[33];
extern const uint64_t G31L[31];

extern const uint64_t T33R[33];
extern const uint64_t T31L[31];

extern const uint64_t N33R[33];
extern const uint64_t N31L[31];

extern const uint64_t* const MS_TAB_33R[ASCII_SIZE];
extern const uint64_t* const MS_TAB_31L[ASCII_SIZE];

extern const uint8_t CONVERT_TAB[ASCII_SIZE];
extern const uint8_t RC_CONVERT_TAB[ASCII_SIZE];

extern const uint64_t DIMER_TAB[4 * 4];
extern const uint64_t TRIMER_TAB[4 * 4 * 4];
extern const uint64_t TETRAMER_TAB[4 * 4 * 4 * 4];

} // namespace btllib

#endif