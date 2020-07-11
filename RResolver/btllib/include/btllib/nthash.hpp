/*
 * nthash.hpp
 * Author: Hamid Mohamadi
 * Genome Sciences Centre,
 * British Columbia Cancer Agency
 */

#ifndef BTLLIB_NTHASH_HPP
#define BTLLIB_NTHASH_HPP

#include <cstdint>
#include <limits>
#include <string>
#include <vector>

namespace btllib {

// offset for the complement base in the random seeds table
const uint8_t cp_off = 0x07;

// shift for gerenerating multiple hash values
const int multishift = 27;

// seed for gerenerating multiple hash values
static const uint64_t multiseed = 0x90b45d39fb6da1fa;

// 64-bit random seeds corresponding to bases and their complements
static const uint64_t seedA = 0x3c8bfbb395c60474;
static const uint64_t seedC = 0x3193c18562a02b4c;
static const uint64_t seedG = 0x20323ed082572324;
static const uint64_t seedT = 0x295549f54be24456;
static const uint64_t seedN = 0x0000000000000000;

static const int ASCII_SIZE = 256;

static const uint64_t seed_tab[ASCII_SIZE] = {
  seedN, seedT, seedN, seedG, seedA, seedA, seedN, seedC, // 0..7
  seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 8..15
  seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 16..23
  seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 24..31
  seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 32..39
  seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 40..47
  seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 48..55
  seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 56..63
  seedN, seedA, seedN, seedC, seedN, seedN, seedN, seedG, // 64..71
  seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 72..79
  seedN, seedN, seedN, seedN, seedT, seedT, seedN, seedN, // 80..87
  seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 88..95
  seedN, seedA, seedN, seedC, seedN, seedN, seedN, seedG, // 96..103
  seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 104..111
  seedN, seedN, seedN, seedN, seedT, seedT, seedN, seedN, // 112..119
  seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 120..127
  seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 128..135
  seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 136..143
  seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 144..151
  seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 152..159
  seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 160..167
  seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 168..175
  seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 176..183
  seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 184..191
  seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 192..199
  seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 200..207
  seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 208..215
  seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 216..223
  seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 224..231
  seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 232..239
  seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 240..247
  seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN  // 248..255
};

static const uint64_t A33r[33] = {
  0x195c60474, 0x12b8c08e9, 0x571811d3,  0xae3023a6,  0x15c60474c, 0xb8c08e99,
  0x171811d32, 0xe3023a65,  0x1c60474ca, 0x18c08e995, 0x11811d32b, 0x3023a657,
  0x60474cae,  0xc08e995c,  0x1811d32b8, 0x1023a6571, 0x474cae3,   0x8e995c6,
  0x11d32b8c,  0x23a65718,  0x474cae30,  0x8e995c60,  0x11d32b8c0, 0x3a657181,
  0x74cae302,  0xe995c604,  0x1d32b8c08, 0x1a6571811, 0x14cae3023, 0x995c6047,
  0x132b8c08e, 0x6571811d,  0xcae3023a
};

static const uint64_t A31l[31] = {
  0x3c8bfbb200000000, 0x7917f76400000000, 0xf22feec800000000,
  0xe45fdd9200000000, 0xc8bfbb2600000000, 0x917f764e00000000,
  0x22feec9e00000000, 0x45fdd93c00000000, 0x8bfbb27800000000,
  0x17f764f200000000, 0x2feec9e400000000, 0x5fdd93c800000000,
  0xbfbb279000000000, 0x7f764f2200000000, 0xfeec9e4400000000,
  0xfdd93c8a00000000, 0xfbb2791600000000, 0xf764f22e00000000,
  0xeec9e45e00000000, 0xdd93c8be00000000, 0xbb27917e00000000,
  0x764f22fe00000000, 0xec9e45fc00000000, 0xd93c8bfa00000000,
  0xb27917f600000000, 0x64f22fee00000000, 0xc9e45fdc00000000,
  0x93c8bfba00000000, 0x27917f7600000000, 0x4f22feec00000000,
  0x9e45fdd800000000
};

static const uint64_t C33r[33] = {
  0x162a02b4c, 0xc5405699,  0x18a80ad32, 0x115015a65, 0x2a02b4cb, 0x54056996,
  0xa80ad32c,  0x15015a658, 0xa02b4cb1,  0x140569962, 0x80ad32c5, 0x1015a658a,
  0x2b4cb15,   0x569962a,   0xad32c54,   0x15a658a8,  0x2b4cb150, 0x569962a0,
  0xad32c540,  0x15a658a80, 0xb4cb1501,  0x169962a02, 0xd32c5405, 0x1a658a80a,
  0x14cb15015, 0x9962a02b,  0x132c54056, 0x658a80ad,  0xcb15015a, 0x1962a02b4,
  0x12c540569, 0x58a80ad3,  0xb15015a6
};

static const uint64_t C31l[31] = {
  0x3193c18400000000, 0x6327830800000000, 0xc64f061000000000,
  0x8c9e0c2200000000, 0x193c184600000000, 0x3278308c00000000,
  0x64f0611800000000, 0xc9e0c23000000000, 0x93c1846200000000,
  0x278308c600000000, 0x4f06118c00000000, 0x9e0c231800000000,
  0x3c18463200000000, 0x78308c6400000000, 0xf06118c800000000,
  0xe0c2319200000000, 0xc184632600000000, 0x8308c64e00000000,
  0x6118c9e00000000,  0xc23193c00000000,  0x1846327800000000,
  0x308c64f000000000, 0x6118c9e000000000, 0xc23193c000000000,
  0x8463278200000000, 0x8c64f0600000000,  0x118c9e0c00000000,
  0x23193c1800000000, 0x4632783000000000, 0x8c64f06000000000,
  0x18c9e0c200000000
};

static const uint64_t G33r[33] = {
  0x82572324,  0x104ae4648, 0x95c8c91,   0x12b91922,  0x25723244,  0x4ae46488,
  0x95c8c910,  0x12b919220, 0x57232441,  0xae464882,  0x15c8c9104, 0xb9192209,
  0x172324412, 0xe4648825,  0x1c8c9104a, 0x191922095, 0x12324412b, 0x46488257,
  0x8c9104ae,  0x11922095c, 0x324412b9,  0x64882572,  0xc9104ae4,  0x1922095c8,
  0x124412b91, 0x48825723,  0x9104ae46,  0x122095c8c, 0x4412b919,  0x88257232,
  0x1104ae464, 0x2095c8c9,  0x412b9192
};

static const uint64_t G31l[31] = {
  0x20323ed000000000, 0x40647da000000000, 0x80c8fb4000000000,
  0x191f68200000000,  0x323ed0400000000,  0x647da0800000000,
  0xc8fb41000000000,  0x191f682000000000, 0x323ed04000000000,
  0x647da08000000000, 0xc8fb410000000000, 0x91f6820200000000,
  0x23ed040600000000, 0x47da080c00000000, 0x8fb4101800000000,
  0x1f68203200000000, 0x3ed0406400000000, 0x7da080c800000000,
  0xfb41019000000000, 0xf682032200000000, 0xed04064600000000,
  0xda080c8e00000000, 0xb410191e00000000, 0x6820323e00000000,
  0xd040647c00000000, 0xa080c8fa00000000, 0x410191f600000000,
  0x820323ec00000000, 0x40647da00000000,  0x80c8fb400000000,
  0x10191f6800000000
};

static const uint64_t T33r[33] = {
  0x14be24456, 0x97c488ad,  0x12f89115a, 0x5f1222b5,  0xbe24456a,  0x17c488ad4,
  0xf89115a9,  0x1f1222b52, 0x1e24456a5, 0x1c488ad4b, 0x189115a97, 0x11222b52f,
  0x24456a5f,  0x488ad4be,  0x9115a97c,  0x1222b52f8, 0x4456a5f1,  0x88ad4be2,
  0x1115a97c4, 0x22b52f89,  0x456a5f12,  0x8ad4be24,  0x115a97c48, 0x2b52f891,
  0x56a5f122,  0xad4be244,  0x15a97c488, 0xb52f8911,  0x16a5f1222, 0xd4be2445,
  0x1a97c488a, 0x152f89115, 0xa5f1222b
};

static const uint64_t T31l[31] = {
  0x295549f400000000, 0x52aa93e800000000, 0xa55527d000000000,
  0x4aaa4fa200000000, 0x95549f4400000000, 0x2aa93e8a00000000,
  0x55527d1400000000, 0xaaa4fa2800000000, 0x5549f45200000000,
  0xaa93e8a400000000, 0x5527d14a00000000, 0xaa4fa29400000000,
  0x549f452a00000000, 0xa93e8a5400000000, 0x527d14aa00000000,
  0xa4fa295400000000, 0x49f452aa00000000, 0x93e8a55400000000,
  0x27d14aaa00000000, 0x4fa2955400000000, 0x9f452aa800000000,
  0x3e8a555200000000, 0x7d14aaa400000000, 0xfa29554800000000,
  0xf452aa9200000000, 0xe8a5552600000000, 0xd14aaa4e00000000,
  0xa295549e00000000, 0x452aa93e00000000, 0x8a55527c00000000,
  0x14aaa4fa00000000
};

static const uint64_t N33r[33] = {
  seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN,
  seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN,
  seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN
};

static const uint64_t N31l[31] = {
  seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN,
  seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN,
  seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN
};

static const uint64_t* ms_tab_33r[ASCII_SIZE] = {
  N33r, T33r, N33r, G33r, A33r, A33r, N33r, C33r, // 0..7
  N33r, N33r, N33r, N33r, N33r, N33r, N33r, N33r, // 8..15
  N33r, N33r, N33r, N33r, N33r, N33r, N33r, N33r, // 16..23
  N33r, N33r, N33r, N33r, N33r, N33r, N33r, N33r, // 24..31
  N33r, N33r, N33r, N33r, N33r, N33r, N33r, N33r, // 32..39
  N33r, N33r, N33r, N33r, N33r, N33r, N33r, N33r, // 40..47
  N33r, N33r, N33r, N33r, N33r, N33r, N33r, N33r, // 48..55
  N33r, N33r, N33r, N33r, N33r, N33r, N33r, N33r, // 56..63
  N33r, A33r, N33r, C33r, N33r, N33r, N33r, G33r, // 64..71
  N33r, N33r, N33r, N33r, N33r, N33r, N33r, N33r, // 72..79
  N33r, N33r, N33r, N33r, T33r, T33r, N33r, N33r, // 80..87
  N33r, N33r, N33r, N33r, N33r, N33r, N33r, N33r, // 88..95
  N33r, A33r, N33r, C33r, N33r, N33r, N33r, G33r, // 96..103
  N33r, N33r, N33r, N33r, N33r, N33r, N33r, N33r, // 104..111
  N33r, N33r, N33r, N33r, T33r, T33r, N33r, N33r, // 112..119
  N33r, N33r, N33r, N33r, N33r, N33r, N33r, N33r, // 120..127
  N33r, N33r, N33r, N33r, N33r, N33r, N33r, N33r, // 128..135
  N33r, N33r, N33r, N33r, N33r, N33r, N33r, N33r, // 136..143
  N33r, N33r, N33r, N33r, N33r, N33r, N33r, N33r, // 144..151
  N33r, N33r, N33r, N33r, N33r, N33r, N33r, N33r, // 152..159
  N33r, N33r, N33r, N33r, N33r, N33r, N33r, N33r, // 160..167
  N33r, N33r, N33r, N33r, N33r, N33r, N33r, N33r, // 168..175
  N33r, N33r, N33r, N33r, N33r, N33r, N33r, N33r, // 176..183
  N33r, N33r, N33r, N33r, N33r, N33r, N33r, N33r, // 184..191
  N33r, N33r, N33r, N33r, N33r, N33r, N33r, N33r, // 192..199
  N33r, N33r, N33r, N33r, N33r, N33r, N33r, N33r, // 200..207
  N33r, N33r, N33r, N33r, N33r, N33r, N33r, N33r, // 208..215
  N33r, N33r, N33r, N33r, N33r, N33r, N33r, N33r, // 216..223
  N33r, N33r, N33r, N33r, N33r, N33r, N33r, N33r, // 224..231
  N33r, N33r, N33r, N33r, N33r, N33r, N33r, N33r, // 232..239
  N33r, N33r, N33r, N33r, N33r, N33r, N33r, N33r, // 240..247
  N33r, N33r, N33r, N33r, N33r, N33r, N33r, N33r  // 248..255
};

static const uint64_t* ms_tab_31l[ASCII_SIZE] = {
  N31l, T31l, N31l, G31l, A31l, A31l, N31l, C31l, // 0..7
  N31l, N31l, N31l, N31l, N31l, N31l, N31l, N31l, // 8..15
  N31l, N31l, N31l, N31l, N31l, N31l, N31l, N31l, // 16..23
  N31l, N31l, N31l, N31l, N31l, N31l, N31l, N31l, // 24..31
  N31l, N31l, N31l, N31l, N31l, N31l, N31l, N31l, // 32..39
  N31l, N31l, N31l, N31l, N31l, N31l, N31l, N31l, // 40..47
  N31l, N31l, N31l, N31l, N31l, N31l, N31l, N31l, // 48..55
  N31l, N31l, N31l, N31l, N31l, N31l, N31l, N31l, // 56..63
  N31l, A31l, N31l, C31l, N31l, N31l, N31l, G31l, // 64..71
  N31l, N31l, N31l, N31l, N31l, N31l, N31l, N31l, // 72..79
  N31l, N31l, N31l, N31l, T31l, T31l, N31l, N31l, // 80..87
  N31l, N31l, N31l, N31l, N31l, N31l, N31l, N31l, // 88..95
  N31l, A31l, N31l, C31l, N31l, N31l, N31l, G31l, // 96..103
  N31l, N31l, N31l, N31l, N31l, N31l, N31l, N31l, // 104..111
  N31l, N31l, N31l, N31l, T31l, T31l, N31l, N31l, // 112..119
  N31l, N31l, N31l, N31l, N31l, N31l, N31l, N31l, // 120..127
  N31l, N31l, N31l, N31l, N31l, N31l, N31l, N31l, // 128..135
  N31l, N31l, N31l, N31l, N31l, N31l, N31l, N31l, // 136..143
  N31l, N31l, N31l, N31l, N31l, N31l, N31l, N31l, // 144..151
  N31l, N31l, N31l, N31l, N31l, N31l, N31l, N31l, // 152..159
  N31l, N31l, N31l, N31l, N31l, N31l, N31l, N31l, // 160..167
  N31l, N31l, N31l, N31l, N31l, N31l, N31l, N31l, // 168..175
  N31l, N31l, N31l, N31l, N31l, N31l, N31l, N31l, // 176..183
  N31l, N31l, N31l, N31l, N31l, N31l, N31l, N31l, // 184..191
  N31l, N31l, N31l, N31l, N31l, N31l, N31l, N31l, // 192..199
  N31l, N31l, N31l, N31l, N31l, N31l, N31l, N31l, // 200..207
  N31l, N31l, N31l, N31l, N31l, N31l, N31l, N31l, // 208..215
  N31l, N31l, N31l, N31l, N31l, N31l, N31l, N31l, // 216..223
  N31l, N31l, N31l, N31l, N31l, N31l, N31l, N31l, // 224..231
  N31l, N31l, N31l, N31l, N31l, N31l, N31l, N31l, // 232..239
  N31l, N31l, N31l, N31l, N31l, N31l, N31l, N31l, // 240..247
  N31l, N31l, N31l, N31l, N31l, N31l, N31l, N31l  // 248..255
};

// rotate "v" to the left 1 position
inline uint64_t
rol1(const uint64_t v)
{
  return (v << 1) | (v >> 63); // NOLINT
}

// rotate "v" to the right by 1 position
inline uint64_t
ror1(const uint64_t v)
{
  return (v >> 1) | (v << 63); // NOLINT
}

// rotate 31-left bits of "v" to the left by "s" positions
inline uint64_t
rol31(const uint64_t v, unsigned s)
{
  s %= 31;                                          // NOLINT
  return ((v << s) | (v >> (31 - s))) & 0x7FFFFFFF; // NOLINT
}

// rotate 33-right bits of "v" to the left by "s" positions
inline uint64_t
rol33(const uint64_t v, unsigned s)
{
  s %= 33;                                           // NOLINT
  return ((v << s) | (v >> (33 - s))) & 0x1FFFFFFFF; // NOLINT
}

// swap bit 0 with bit 33 in "v"
inline uint64_t
swapbits033(const uint64_t v)
{
  uint64_t x = (v ^ (v >> 33)) & 1; // NOLINT
  return v ^ (x | (x << 33));       // NOLINT
}

// swap bit 32 with bit 63 in "v"
inline uint64_t
swapbits3263(const uint64_t v)
{
  uint64_t x = ((v >> 32) ^ (v >> 63)) & 1; // NOLINT
  return v ^ ((x << 32) | (x << 63));       // NOLINT
}

// forward-strand hash value of the base kmer, i.e. fhval(kmer_0)
inline uint64_t
NTF64(const char* kmer_seq, const unsigned k)
{
  uint64_t h_val = 0;
  for (unsigned i = 0; i < k; i++) {
    h_val = rol1(h_val);
    h_val = swapbits033(h_val);
    h_val ^= seed_tab[(unsigned char)kmer_seq[i]];
  }
  return h_val;
}

// reverse-strand hash value of the base kmer, i.e. rhval(kmer_0)
inline uint64_t
NTR64(const char* kmer_seq, const unsigned k)
{
  uint64_t h_val = 0;
  for (unsigned i = 0; i < k; i++) {
    h_val = rol1(h_val);
    h_val = swapbits033(h_val);
    h_val ^= seed_tab[(unsigned char)kmer_seq[k - 1 - i] & cp_off];
  }
  return h_val;
}

// forward-strand ntHash for sliding k-mers
inline uint64_t
NTF64(const uint64_t fh_val,
      const unsigned k,
      const unsigned char char_out,
      const unsigned char char_in)
{
  uint64_t h_val = rol1(fh_val);
  h_val = swapbits033(h_val);
  h_val ^= seed_tab[char_in];
  h_val ^=
    (ms_tab_31l[char_out][k % 31] | ms_tab_33r[char_out][k % 33]); // NOLINT
  return h_val;
}

// reverse-complement ntHash for sliding k-mers
inline uint64_t
NTR64(const uint64_t rh_val,
      const unsigned k,
      const unsigned char char_out,
      const unsigned char char_in)
{
  uint64_t h_val = rh_val ^ (ms_tab_31l[char_in & cp_off][k % 31] | // NOLINT
                             ms_tab_33r[char_in & cp_off][k % 33]); // NOLINT
  h_val ^= seed_tab[char_out & cp_off];
  h_val = ror1(h_val);
  h_val = swapbits3263(h_val);
  return h_val;
}

// canonical ntBase
inline uint64_t
NTC64(const char* kmer_seq, const unsigned k)
{
  uint64_t fh_val = 0, rh_val = 0;
  fh_val = NTF64(kmer_seq, k);
  rh_val = NTR64(kmer_seq, k);
  return (rh_val < fh_val) ? rh_val : fh_val;
}

// canonical ntHash
inline uint64_t
NTC64(const char* kmer_seq,
      const unsigned k,
      uint64_t& fh_val,
      uint64_t& rh_val)
{
  fh_val = NTF64(kmer_seq, k);
  rh_val = NTR64(kmer_seq, k);
  return (rh_val < fh_val) ? rh_val : fh_val;
}

// canonical ntHash for sliding k-mers
inline uint64_t
NTC64(const unsigned char char_out,
      const unsigned char char_in,
      const unsigned k,
      uint64_t& fh_val,
      uint64_t& rh_val)
{
  fh_val = NTF64(fh_val, k, char_out, char_in);
  rh_val = NTR64(rh_val, k, char_out, char_in);
  return (rh_val < fh_val) ? rh_val : fh_val;
}

// forward-strand ntHash for sliding k-mers to the left
inline uint64_t
NTF64L(const uint64_t rh_val,
       const unsigned k,
       const unsigned char char_out,
       const unsigned char char_in)
{
  uint64_t h_val = rh_val ^ (ms_tab_31l[char_in][k % 31] | // NOLINT
                             ms_tab_33r[char_in][k % 33]); // NOLINT
  h_val ^= seed_tab[char_out];
  h_val = ror1(h_val);
  h_val = swapbits3263(h_val);
  return h_val;
}

// reverse-complement ntHash for sliding k-mers to the left
inline uint64_t
NTR64L(const uint64_t fh_val,
       const unsigned k,
       const unsigned char char_out,
       const unsigned char char_in)
{
  uint64_t h_val = rol1(fh_val);
  h_val = swapbits033(h_val);
  h_val ^= seed_tab[char_in & cp_off];
  h_val ^= (ms_tab_31l[char_out & cp_off][k % 31] | // NOLINT
            ms_tab_33r[char_out & cp_off][k % 33]); // NOLINT
  return h_val;
}

// canonical ntHash for sliding k-mers to the left
inline uint64_t
NTC64L(const unsigned char char_out,
       const unsigned char char_in,
       const unsigned k,
       uint64_t& fh_val,
       uint64_t& rh_val)
{
  fh_val = NTF64L(fh_val, k, char_out, char_in);
  rh_val = NTR64L(rh_val, k, char_out, char_in);
  return (rh_val < fh_val) ? rh_val : fh_val;
}

// ntBase with seeding option
inline uint64_t
NTF64(const char* kmer_seq, const unsigned k, const unsigned seed)
{
  uint64_t h_val = NTF64(kmer_seq, k);
  if (seed == 0) {
    return h_val;
  }
  h_val *= seed ^ k * multiseed;
  h_val ^= h_val >> multishift;
  return h_val;
}

// canonical ntBase with seeding option
inline uint64_t
NTC64(const char* kmer_seq, const unsigned k, const unsigned seed)
{
  uint64_t h_val = NTC64(kmer_seq, k);
  if (seed == 0) {
    return h_val;
  }
  h_val *= seed ^ k * multiseed;
  h_val ^= h_val >> multishift;
  return h_val;
}

// multihash ntHash, ntBase
inline void
NTM64(const char* kmer_seq, const unsigned k, const unsigned m, uint64_t* h_val)
{
  uint64_t b_val = 0, t_val = 0;
  b_val = NTF64(kmer_seq, k);
  h_val[0] = b_val;
  for (unsigned i = 1; i < m; i++) {
    t_val = b_val * (i ^ k * multiseed);
    t_val ^= t_val >> multishift;
    h_val[i] = t_val;
  }
}

// one extra hash for given base hash
inline uint64_t
NTE64(const uint64_t h_val, const unsigned k, const unsigned i)
{
  uint64_t t_val = h_val;
  t_val *= (i ^ k * multiseed);
  t_val ^= t_val >> multishift;
  return t_val;
}

// multihash ntHash for sliding k-mers
inline void
NTM64(const unsigned char char_out,
      const unsigned char char_in,
      const unsigned k,
      const unsigned m,
      uint64_t* h_val)
{
  uint64_t b_val = 0, t_val = 0;
  b_val = NTF64(h_val[0], k, char_out, char_in);
  h_val[0] = b_val;
  for (unsigned i = 1; i < m; i++) {
    t_val = b_val * (i ^ k * multiseed);
    t_val ^= t_val >> multishift;
    h_val[i] = t_val;
  }
}

// canonical multihash ntBase
inline void
NTMC64(const char* kmer_seq,
       const unsigned k,
       const unsigned m,
       uint64_t* h_val)
{
  uint64_t b_val = 0, t_val = 0;
  b_val = NTC64(kmer_seq, k);
  h_val[0] = b_val;
  for (unsigned i = 1; i < m; i++) {
    t_val = b_val * (i ^ k * multiseed);
    t_val ^= t_val >> multishift;
    h_val[i] = t_val;
  }
}

// canonical multihash ntHash
inline void
NTMC64(const char* kmer_seq,
       const unsigned k,
       const unsigned m,
       uint64_t& fh_val,
       uint64_t& rh_val,
       uint64_t* h_val)
{
  uint64_t b_val = 0, t_val = 0;
  b_val = NTC64(kmer_seq, k, fh_val, rh_val);
  h_val[0] = b_val;
  for (unsigned i = 1; i < m; i++) {
    t_val = b_val * (i ^ k * multiseed);
    t_val ^= t_val >> multishift;
    h_val[i] = t_val;
  }
}

// canonical multihash ntHash for sliding k-mers
inline void
NTMC64(const unsigned char char_out,
       const unsigned char char_in,
       const unsigned k,
       const unsigned m,
       uint64_t& fh_val,
       uint64_t& rh_val,
       uint64_t* h_val)
{
  uint64_t b_val = 0, t_val = 0;
  b_val = NTC64(char_out, char_in, k, fh_val, rh_val);
  h_val[0] = b_val;
  for (unsigned i = 1; i < m; i++) {
    t_val = b_val * (i ^ k * multiseed);
    t_val ^= t_val >> multishift;
    h_val[i] = t_val;
  }
}

/*
 * ignoring k-mers containing nonACGT using ntHash function
 */

// canonical ntBase
inline bool
NTC64(const char* kmer_seq, const unsigned k, uint64_t& h_val, unsigned& locN)
{
  h_val = 0;
  locN = 0;
  uint64_t fh_val = 0, rh_val = 0;
  for (int i = int(k - 1); i >= 0; i--) {
    if (seed_tab[(unsigned char)kmer_seq[i]] == seedN) {
      locN = i;
      return false;
    }
    fh_val = rol1(fh_val);
    fh_val = swapbits033(fh_val);
    fh_val ^= seed_tab[(unsigned char)kmer_seq[k - 1 - i]];

    rh_val = rol1(rh_val);
    rh_val = swapbits033(rh_val);
    rh_val ^= seed_tab[(unsigned char)kmer_seq[i] & cp_off];
  }
  h_val = (rh_val < fh_val) ? rh_val : fh_val;
  return true;
}

// canonical multihash ntBase
inline bool
NTMC64(const char* kmer_seq,
       const unsigned k,
       const unsigned m,
       unsigned& locN,
       uint64_t* h_val)
{
  uint64_t b_val = 0, t_val = 0, fh_val = 0, rh_val = 0;
  locN = 0;
  for (int i = int(k - 1); i >= 0; i--) {
    if (seed_tab[(unsigned char)kmer_seq[i]] == seedN) {
      locN = i;
      return false;
    }
    fh_val = rol1(fh_val);
    fh_val = swapbits033(fh_val);
    fh_val ^= seed_tab[(unsigned char)kmer_seq[k - 1 - i]];

    rh_val = rol1(rh_val);
    rh_val = swapbits033(rh_val);
    rh_val ^= seed_tab[(unsigned char)kmer_seq[i] & cp_off];
  }
  b_val = (rh_val < fh_val) ? rh_val : fh_val;
  h_val[0] = b_val;
  for (unsigned i = 1; i < m; i++) {
    t_val = b_val * (i ^ k * multiseed);
    t_val ^= t_val >> multishift;
    h_val[i] = t_val;
  }
  return true;
}

// canonical ntHash
inline bool
NTC64(const char* kmer_seq,
      const unsigned k,
      uint64_t& fh_val,
      uint64_t& rh_val,
      uint64_t& h_val,
      unsigned& locN)
{
  h_val = fh_val = rh_val = 0;
  locN = 0;
  for (int i = int(k - 1); i >= 0; i--) {
    if (seed_tab[(unsigned char)kmer_seq[i]] == seedN) {
      locN = i;
      return false;
    }
    fh_val = rol1(fh_val);
    fh_val = swapbits033(fh_val);
    fh_val ^= seed_tab[(unsigned char)kmer_seq[k - 1 - i]];

    rh_val = rol1(rh_val);
    rh_val = swapbits033(rh_val);
    rh_val ^= seed_tab[(unsigned char)kmer_seq[i] & cp_off];
  }
  h_val = (rh_val < fh_val) ? rh_val : fh_val;
  return true;
}

// canonical multihash ntHash
inline bool
NTMC64(const char* kmer_seq,
       const unsigned k,
       const unsigned m,
       uint64_t& fh_val,
       uint64_t& rh_val,
       unsigned& locN,
       uint64_t* h_val)
{
  fh_val = rh_val = 0;
  uint64_t b_val = 0, t_val = 0;
  locN = 0;
  for (int i = int(k - 1); i >= 0; i--) {
    if (seed_tab[(unsigned char)kmer_seq[i]] == seedN) {
      locN = i;
      return false;
    }
    fh_val = rol1(fh_val);
    fh_val = swapbits033(fh_val);
    fh_val ^= seed_tab[(unsigned char)kmer_seq[k - 1 - i]];

    rh_val = rol1(rh_val);
    rh_val = swapbits033(rh_val);
    rh_val ^= seed_tab[(unsigned char)kmer_seq[i] & cp_off];
  }
  b_val = (rh_val < fh_val) ? rh_val : fh_val;
  h_val[0] = b_val;
  for (unsigned i = 1; i < m; i++) {
    t_val = b_val * (i ^ k * multiseed);
    t_val ^= t_val >> multishift;
    h_val[i] = t_val;
  }
  return true;
}

// strand-aware canonical multihash ntHash
inline bool
NTMC64(const char* kmer_seq,
       const unsigned k,
       const unsigned m,
       uint64_t& fh_val,
       uint64_t& rh_val,
       unsigned& locN,
       uint64_t* h_val,
       bool& hStn)
{
  fh_val = rh_val = 0;
  uint64_t b_val = 0, t_val = 0;
  locN = 0;
  for (int i = int(k - 1); i >= 0; i--) {
    if (seed_tab[(unsigned char)kmer_seq[i]] == seedN) {
      locN = i;
      return false;
    }
    fh_val = rol1(fh_val);
    fh_val = swapbits033(fh_val);
    fh_val ^= seed_tab[(unsigned char)kmer_seq[k - 1 - i]];

    rh_val = rol1(rh_val);
    rh_val = swapbits033(rh_val);
    rh_val ^= seed_tab[(unsigned char)kmer_seq[i] & cp_off];
  }
  hStn = rh_val < fh_val;
  b_val = hStn ? rh_val : fh_val;
  h_val[0] = b_val;
  for (unsigned i = 1; i < m; i++) {
    t_val = b_val * (i ^ k * multiseed);
    t_val ^= t_val >> multishift;
    h_val[i] = t_val;
  }
  return true;
}

// starnd-aware canonical multihash ntHash for sliding k-mers
inline void
NTMC64(const unsigned char char_out,
       const unsigned char char_in,
       const unsigned k,
       const unsigned m,
       uint64_t& fh_val,
       uint64_t& rh_val,
       uint64_t* h_val,
       bool& hStn)
{
  uint64_t b_val = 0, t_val = 0;
  b_val = NTC64(char_out, char_in, k, fh_val, rh_val);
  hStn = rh_val < fh_val;
  h_val[0] = b_val;
  for (unsigned i = 1; i < m; i++) {
    t_val = b_val * (i ^ k * multiseed);
    t_val ^= t_val >> multishift;
    h_val[i] = t_val;
  }
}

// masking canonical ntHash using spaced seed pattern
inline uint64_t
maskHash(uint64_t& fk_val,
         uint64_t& rk_val,
         const char* seed_seq,
         const char* kmer_seq,
         const unsigned k)
{
  uint64_t fs_val = fk_val, rs_val = rk_val;
  for (unsigned i = 0; i < k; i++) {
    if (seed_seq[i] != '1') {
      fs_val ^=
        (ms_tab_31l[(unsigned char)kmer_seq[i]][(k - 1 - i) % 31] | // NOLINT
         ms_tab_33r[(unsigned char)kmer_seq[i]][(k - 1 - i) % 33]); // NOLINT
      rs_val ^=
        (ms_tab_31l[(unsigned char)kmer_seq[i] & cp_off][i % 31] | // NOLINT
         ms_tab_33r[(unsigned char)kmer_seq[i] & cp_off][i % 33]); // NOLINT
    }
  }
  return (rs_val < fs_val) ? rs_val : fs_val;
}

// spaced seed ntHash for base kmer, i.e. fhval(kmer_0)
inline uint64_t
NTS64(const char* kmer_seq,
      const std::vector<bool>& seed,
      const unsigned k,
      uint64_t& h_val)
{
  h_val = 0;
  uint64_t sVal = 0;
  for (unsigned i = 0; i < k; i++) {
    h_val = rol1(h_val);
    h_val = swapbits033(h_val);
    sVal = h_val;
    h_val ^= seed_tab[(unsigned char)kmer_seq[i]];
    if (seed[i]) {
      sVal = h_val;
    }
  }
  return sVal;
}

// spaced seed ntHash for sliding k-mers
inline uint64_t
NTS64(const char* kmer_seq,
      const std::vector<bool>& seed,
      const unsigned char char_out,
      const unsigned char char_in,
      const unsigned k,
      uint64_t& h_val)
{
  h_val = NTF64(h_val, k, char_out, char_in);
  uint64_t sVal = h_val;
  for (unsigned i = 0; i < k; i++) {
    if (!seed[i]) {
      sVal ^= (ms_tab_31l[(unsigned char)kmer_seq[i]][k % 31] | // NOLINT
               ms_tab_33r[(unsigned char)kmer_seq[i]][k % 33]); // NOLINT
    }
  }
  return sVal;
}

// strand-aware multihash spaced seed ntHash
inline bool
NTMS64(const char* kmer_seq,
       const std::vector<std::vector<unsigned>>& seed_seq,
       const unsigned k,
       const unsigned m,
       uint64_t& fh_val,
       uint64_t& rh_val,
       unsigned& locN,
       uint64_t* h_val,
       bool* hStn)
{
  fh_val = rh_val = 0;
  locN = 0;
  for (int i = int(k - 1); i >= 0; i--) {
    if (seed_tab[(unsigned char)kmer_seq[i]] == seedN) {
      locN = i;
      return false;
    }
    fh_val = rol1(fh_val);
    fh_val = swapbits033(fh_val);
    fh_val ^= seed_tab[(unsigned char)kmer_seq[k - 1 - i]];

    rh_val = rol1(rh_val);
    rh_val = swapbits033(rh_val);
    rh_val ^= seed_tab[(unsigned char)kmer_seq[i] & cp_off];
  }

  for (unsigned j = 0; j < m; j++) {
    uint64_t fs_val = fh_val, rs_val = rh_val;
    for (const auto& seed_pos : seed_seq[j]) {
      fs_val ^= (ms_tab_31l[(unsigned char)kmer_seq[seed_pos]]
                           [(k - 1 - seed_pos) % 31] | // NOLINT
                 ms_tab_33r[(unsigned char)kmer_seq[seed_pos]]
                           [(k - 1 - seed_pos) % 33]); // NOLINT
      rs_val ^= (ms_tab_31l[(unsigned char)kmer_seq[seed_pos] & cp_off]
                           [seed_pos % 31] | // NOLINT
                 ms_tab_33r[(unsigned char)kmer_seq[seed_pos] & cp_off]
                           [seed_pos % 33]); // NOLINT
    }
    hStn[j] = rs_val < fs_val;
    h_val[j] = hStn[j] ? rs_val : fs_val;
  }
  return true;
}

// strand-aware multihash spaced seed ntHash for sliding k-mers
inline void
NTMS64(const char* kmer_seq,
       const std::vector<std::vector<unsigned>>& seed_seq,
       const unsigned char char_out,
       const unsigned char char_in,
       const unsigned k,
       const unsigned m,
       uint64_t& fh_val,
       uint64_t& rh_val,
       uint64_t* h_val,
       bool* hStn)
{
  fh_val = NTF64(fh_val, k, char_out, char_in);
  rh_val = NTR64(rh_val, k, char_out, char_in);
  for (unsigned j = 0; j < m; j++) {
    uint64_t fs_val = fh_val, rs_val = rh_val;
    for (const auto& seed_pos : seed_seq[j]) {
      fs_val ^= (ms_tab_31l[(unsigned char)kmer_seq[seed_pos]]
                           [(k - 1 - seed_pos) % 31] | // NOLINT
                 ms_tab_33r[(unsigned char)kmer_seq[seed_pos]]
                           [(k - 1 - seed_pos) % 33]); // NOLINT
      rs_val ^= (ms_tab_31l[(unsigned char)kmer_seq[seed_pos] & cp_off]
                           [seed_pos % 31] | // NOLINT
                 ms_tab_33r[(unsigned char)kmer_seq[seed_pos] & cp_off]
                           [seed_pos % 33]); // NOLINT
      ;
    }
    hStn[j] = rs_val < fs_val;
    h_val[j] = hStn[j] ? rs_val : fs_val;
  }
}

// Multi spaced seed ntHash with multiple hashes per seed
inline bool
NTMSM64(const char* kmer_seq,
        const std::vector<std::vector<unsigned>>& seed_seq,
        const unsigned k,
        const unsigned m,
        const unsigned m2,
        uint64_t& fh_val,
        uint64_t& rh_val,
        unsigned& locN,
        uint64_t* h_val)
{
  fh_val = rh_val = 0;
  locN = 0;
  for (int i = int(k - 1); i >= 0; i--) {
    if (seed_tab[(unsigned char)kmer_seq[i]] == seedN) {
      locN = i;
      return false;
    }
    fh_val = rol1(fh_val);
    fh_val = swapbits033(fh_val);
    fh_val ^= seed_tab[(unsigned char)kmer_seq[k - 1 - i]];

    rh_val = rol1(rh_val);
    rh_val = swapbits033(rh_val);
    rh_val ^= seed_tab[(unsigned char)kmer_seq[i] & cp_off];
  }

  for (unsigned j = 0; j < m; j++) {
    uint64_t fs_val = fh_val, rs_val = rh_val;
    for (const auto& seed_pos : seed_seq[j]) {
      fs_val ^= (ms_tab_31l[(unsigned char)kmer_seq[seed_pos]]
                           [(k - 1 - seed_pos) % 31] | // NOLINT
                 ms_tab_33r[(unsigned char)kmer_seq[seed_pos]]
                           [(k - 1 - seed_pos) % 33]); // NOLINT
      rs_val ^= (ms_tab_31l[(unsigned char)kmer_seq[seed_pos] & cp_off]
                           [seed_pos % 31] | // NOLINT
                 ms_tab_33r[(unsigned char)kmer_seq[seed_pos] & cp_off]
                           [seed_pos % 33]); // NOLINT
    }
    h_val[j * m2] = rs_val < fs_val ? rs_val : fs_val;
    for (unsigned j2 = 1; j2 < m2; j2++) {
      uint64_t t_val = h_val[j * m2] * (j2 ^ k * multiseed);
      t_val ^= t_val >> multishift;
      h_val[j * m2 + j2] = t_val;
    }
  }
  return true;
}

// Multi spaced seed ntHash for sliding k-mers with multiple hashes per seed
inline void
NTMSM64(const char* kmer_seq,
        const std::vector<std::vector<unsigned>>& seed_seq,
        const unsigned char char_out,
        const unsigned char char_in,
        const unsigned k,
        const unsigned m,
        const unsigned m2,
        uint64_t& fh_val,
        uint64_t& rh_val,
        uint64_t* h_val)
{
  fh_val = NTF64(fh_val, k, char_out, char_in);
  rh_val = NTR64(rh_val, k, char_out, char_in);
  for (unsigned j = 0; j < m; j++) {
    uint64_t fs_val = fh_val, rs_val = rh_val;
    for (const auto& seed_pos : seed_seq[j]) {
      fs_val ^= (ms_tab_31l[(unsigned char)kmer_seq[seed_pos]]
                           [(k - 1 - seed_pos) % 31] | // NOLINT
                 ms_tab_33r[(unsigned char)kmer_seq[seed_pos]]
                           [(k - 1 - seed_pos) % 33]); // NOLINT
      rs_val ^= (ms_tab_31l[(unsigned char)kmer_seq[seed_pos] & cp_off]
                           [seed_pos % 31] | // NOLINT
                 ms_tab_33r[(unsigned char)kmer_seq[seed_pos] & cp_off]
                           [seed_pos % 33]); // NOLINT
    }
    h_val[j * m2] = rs_val < fs_val ? rs_val : fs_val;
    for (unsigned j2 = 1; j2 < m2; j2++) {
      uint64_t t_val = h_val[j * m2] * (j2 ^ k * multiseed);
      t_val ^= t_val >> multishift;
      h_val[j * m2 + j2] = t_val;
    }
  }
}

class NtHash;
class SeedNtHash;
using SpacedSeed = std::vector<unsigned>;
static std::vector<SpacedSeed>
parse_seeds(const std::vector<std::string>& seed_strings);

/**
 * Iterate over hash values for k-mers in a
 * given DNA sequence.
 *
 * This implementation uses ntHash
 * function to efficiently calculate
 * hash values for successive k-mers.
 */
class NtHash
{

public:
  /**
   * Constructor.
   * @param seq DNA sequence to be hashed
   * @param seq_len length of seq
   * @param k k-mer size
   * @param hash_num number of hashes
   */
  NtHash(const char* seq, size_t seq_len, unsigned k, unsigned hash_num);

  /**
   * Constructor.
   * @param seq DNA sequence to be hashed
   * @param k k-mer size
   * @param hash_num number of hashes
   */
  NtHash(const std::string& seq, unsigned k, unsigned hash_num);

  /**
   * Calculate the next hash value
   * @return true on success and false otherwise
   */
  bool roll();

  const uint64_t* hashes() const;

  size_t get_pos() const { return pos; }
  unsigned get_k() const { return k; }
  unsigned get_hash_num() const { return hash_num; }

protected:
  /** Initialize internal state of iterator */
  bool init();

  const char* seq;
  const size_t seq_len;
  const unsigned k;
  const unsigned hash_num;
  size_t pos = 0;
  std::vector<uint64_t> hashes_vector;
  uint64_t forward_hash = 0;
  uint64_t reverse_hash = 0;
};

class SeedNtHash : public NtHash
{

public:
  SeedNtHash(const char* seq,
             size_t seq_len,
             unsigned k,
             const std::vector<SpacedSeed>& seeds,
             unsigned hash_num_per_seed);
  SeedNtHash(const std::string& seq,
             unsigned k,
             const std::vector<SpacedSeed>& seeds,
             unsigned hash_num_per_seed);
  SeedNtHash(const char* seq,
             size_t seq_len,
             unsigned k,
             const std::vector<std::string>& seeds,
             unsigned hash_num_per_seed);
  SeedNtHash(const std::string& seq,
             unsigned k,
             const std::vector<std::string>& seeds,
             unsigned hash_num_per_seed);

  unsigned get_hash_num_per_seed() const { return hash_num_per_seed; }

  bool roll();

private:
  bool init();

  const unsigned hash_num_per_seed;
  std::vector<SpacedSeed> seeds;
};

inline NtHash::NtHash(const char* seq,
                      size_t seq_len,
                      unsigned k,
                      unsigned hash_num)
  : seq(seq)
  , seq_len(seq_len)
  , k(k)
  , hash_num(hash_num)
{
  hashes_vector.resize(hash_num);
}

inline NtHash::NtHash(const std::string& seq, unsigned k, unsigned hash_num)
  : NtHash(seq.c_str(), seq.size(), k, hash_num)
{}

inline SeedNtHash::SeedNtHash(const char* seq,
                              size_t seq_len,
                              unsigned k,
                              const std::vector<SpacedSeed>& seeds,
                              unsigned hash_num_per_seed)
  : NtHash(seq, seq_len, k, seeds.size() * hash_num_per_seed)
  , hash_num_per_seed(hash_num_per_seed)
  , seeds(seeds)
{}

inline SeedNtHash::SeedNtHash(const std::string& seq,
                              unsigned k,
                              const std::vector<SpacedSeed>& seeds,
                              unsigned hash_num_per_seed)
  : NtHash(seq, k, seeds.size() * hash_num_per_seed)
  , hash_num_per_seed(hash_num_per_seed)
  , seeds(seeds)
{}

inline SeedNtHash::SeedNtHash(const char* seq,
                              size_t seq_len,
                              unsigned k,
                              const std::vector<std::string>& seeds,
                              unsigned hash_num_per_seed)
  : NtHash(seq, seq_len, k, seeds.size() * hash_num_per_seed)
  , hash_num_per_seed(hash_num_per_seed)
  , seeds(parse_seeds(seeds))
{}

inline SeedNtHash::SeedNtHash(const std::string& seq,
                              unsigned k,
                              const std::vector<std::string>& seeds,
                              unsigned hash_num_per_seed)
  : NtHash(seq, k, seeds.size() * hash_num_per_seed)
  , hash_num_per_seed(hash_num_per_seed)
  , seeds(parse_seeds(seeds))
{}

static std::vector<SpacedSeed>
parse_seeds(const std::vector<std::string>& seed_strings)
{
  std::vector<SpacedSeed> seed_set;
  for (const auto& seed_string : seed_strings) {
    SpacedSeed seed;
    size_t pos = 0;
    for (const auto& c : seed_string) {
      if (c != '1') {
        seed.push_back(pos);
      }
      ++pos;
    }
    seed_set.push_back(seed);
  }
  return seed_set;
}

// NOLINTNEXTLINE
#define NT_HASH_INIT(CLASS, NTHASH_CALL)                                       \
  inline bool CLASS::init()                                                    \
  {                                                                            \
    if (k > seq_len) {                                                         \
      pos = std::numeric_limits<std::size_t>::max();                           \
      return false;                                                            \
    }                                                                          \
    unsigned posN = 0;                                                         \
    while ((pos < seq_len - k + 1) && !(NTHASH_CALL)) {                        \
      pos += posN + 1;                                                         \
    }                                                                          \
    if (pos > seq_len - k) {                                                   \
      pos = std::numeric_limits<std::size_t>::max();                           \
      return false;                                                            \
    }                                                                          \
    ++pos;                                                                     \
    return true;                                                               \
  }

// NOLINTNEXTLINE
#define NT_HASH_ROLL(CLASS, NTHASH_CALL)                                       \
  inline bool CLASS::roll()                                                    \
  {                                                                            \
    if (pos == 0) {                                                            \
      return init();                                                           \
    }                                                                          \
    if (pos > seq_len - k) {                                                   \
      return false;                                                            \
    }                                                                          \
    if (seed_tab[(unsigned char)(seq[pos + k - 1])] == seedN) {                \
      pos += k;                                                                \
      return init();                                                           \
    }                                                                          \
    (NTHASH_CALL);                                                             \
    ++pos;                                                                     \
    return true;                                                               \
  }

NT_HASH_INIT(NtHash,
             NTMC64(seq + pos,
                    k,
                    hash_num,
                    forward_hash,
                    reverse_hash,
                    posN,
                    hashes_vector.data()))
NT_HASH_ROLL(NtHash,
             NTMC64(seq[pos - 1],
                    seq[pos - 1 + k],
                    k,
                    hash_num,
                    forward_hash,
                    reverse_hash,
                    hashes_vector.data()))

NT_HASH_INIT(SeedNtHash,
             NTMSM64(seq + pos,
                     seeds,
                     k,
                     seeds.size(),
                     hash_num_per_seed,
                     forward_hash,
                     reverse_hash,
                     posN,
                     hashes_vector.data()))
NT_HASH_ROLL(SeedNtHash,
             NTMSM64(seq + pos,
                     seeds,
                     seq[pos - 1],
                     seq[pos - 1 + k],
                     k,
                     seeds.size(),
                     hash_num_per_seed,
                     forward_hash,
                     reverse_hash,
                     hashes_vector.data()))

#undef NT_HASH_INIT
#undef NT_HASH_ROLL

inline const uint64_t*
NtHash::hashes() const
{
  return hashes_vector.data();
}

} // namespace btllib

#endif
