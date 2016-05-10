#ifndef ROLLING_HASH_H
#define ROLLING_HASH_H

#include <stdint.h>

// offset for the complement base in the random seeds table
const int cpOff = -20;

// shift for gerenerating multiple hash values
const int varShift = 27;

// seed for gerenerating multiple hash values
const uint64_t varSeed = 0x90b45d39fb6da1fa;

// 64-bit random seed for each base
static const uint64_t seedN = 0x0000000000000000;
static const uint64_t seedA = 0x3c8bfbb395c60474;
static const uint64_t seedC = 0x3193c18562a02b4c;
static const uint64_t seedG = 0x20323ed082572324;
static const uint64_t seedT = 0x295549f54be24456;

// 64-bit random seed table corresponding to bases and their complements
static const uint64_t seedTab[256] = {
    seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 0..7
    seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 8..15
    seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 16..23
    seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 24..31
    seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 32..39
    seedN, seedN, seedN, seedN, seedN, seedT, seedN, seedG, // 40..47
    seedN, seedN, seedN, seedC, seedN, seedN, seedN, seedN, // 48..55
    seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 56..63
    seedA, seedA, seedN, seedC, seedN, seedN, seedN, seedG, // 64..71
    seedN, seedN, seedN, seedN, seedN, seedT, seedN, seedG, // 72..79
    seedN, seedN, seedN, seedC, seedT, seedN, seedN, seedN, // 80..87
    seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 88..95
    seedA, seedA, seedN, seedC, seedN, seedN, seedN, seedG, // 96..103
    seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 104..111
    seedN, seedN, seedN, seedN, seedT, seedN, seedN, seedN, // 112..119
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

// rotate "v" to the left by "s" positions
inline uint64_t rol(const uint64_t v, const int s) {
    return (v << s) | (v >> (64 - s));
}

// rotate "v" to the right by "s" positions
inline uint64_t ror(const uint64_t v, const int s) {
    return (v >> s) | (v << (64 - s));
}

// forward-strand hash value of the base kmer, i.e. fhval(kmer_0)
inline uint64_t getFhval(const char * kmerSeq, const unsigned k) {
    uint64_t hVal=0;
    for(unsigned i=0; i<k; i++)
        hVal ^= rol(seedTab[(unsigned char)kmerSeq[i]], k-1-i);
    return hVal;
}

// reverse-strand hash value of the base kmer, i.e. rhval(kmer_0)
inline uint64_t getRhval(const char * kmerSeq, const unsigned k) {
    uint64_t hVal=0;
    for(unsigned i=0; i<k; i++)
        hVal ^= rol(seedTab[(unsigned char)kmerSeq[i]+cpOff], i);
    return hVal;
}

// cannonical hash value of the base kmer, i.e. rhval(kmer_0)
inline uint64_t getChval(const char * kmerSeq, const unsigned k) {
    uint64_t fhVal = getFhval(kmerSeq, k);
    uint64_t rhVal = getRhval(kmerSeq, k);
    return rhVal^fhVal;
}

// initialize forward-strand hash value of the first kmer, i.e. fhval(kmer_0)
inline uint64_t initHashes(const char * kmerSeq, const unsigned k) {
    return getFhval(kmerSeq, k);
}

// initialize cannonical hash value of the first kmer, i.e. chval(kmer_0)
inline uint64_t initHashes(const char * kmerSeq, const unsigned k, uint64_t& fhVal, uint64_t& rhVal) {
    fhVal = getFhval(kmerSeq, k);
    rhVal = getRhval(kmerSeq, k);
    return rhVal^fhVal;
}

// recursive forward-strand hash value for next k-mer
inline uint64_t rollHashesRight(const uint64_t fhVal, const unsigned char charOut, const unsigned char charIn, const unsigned k) {
    return(rol(fhVal, 1) ^ rol(seedTab[charOut], k) ^ seedTab[charIn]);
}

// recursive cannonical hash value for next k-mer
inline uint64_t rollHashesRight(uint64_t& fhVal, uint64_t& rhVal, const unsigned char charOut, const unsigned char charIn, const unsigned k) {
    fhVal = rol(fhVal, 1) ^ rol(seedTab[charOut], k) ^ seedTab[charIn];
    rhVal = ror(rhVal, 1) ^ ror(seedTab[charOut+cpOff], 1) ^ rol(seedTab[charIn+cpOff], k-1);
    return rhVal^fhVal;
}

// recursive forward-strand hash value for prev k-mer
inline uint64_t rollHashesLeft(const uint64_t fhVal, const unsigned char charIn, const unsigned char charOut, const unsigned k) {
    return(ror(fhVal, 1) ^ ror(seedTab[charOut], 1) ^ rol(seedTab[charIn], k-1));
}

// recursive canonical hash value for prev k-mer
inline uint64_t rollHashesLeft(uint64_t& fhVal, uint64_t& rhVal, const unsigned char charIn, const unsigned char charOut, const unsigned k) {
    fhVal = ror(fhVal, 1) ^ ror(seedTab[charOut], 1) ^ rol(seedTab[charIn], k-1);
    rhVal = rol(rhVal, 1) ^ rol(seedTab[charOut+cpOff], k) ^ seedTab[charIn+cpOff];
    return rhVal^fhVal;
}

// change a single base and update forward-strand hash value accordingly
inline uint64_t setBase(uint64_t fhVal, char* kmerSeq, unsigned pos, char base, unsigned k)
{
    fhVal ^= rol(seedTab[(unsigned char)kmerSeq[pos]], k-1-pos);
    kmerSeq[pos] = base;
    fhVal ^= rol(seedTab[(unsigned char)kmerSeq[pos]], k-1-pos);
    return fhVal;
}

// change a single base and update hash values accordingly
inline uint64_t setBase(uint64_t& fhVal, uint64_t& rhVal, char* kmerSeq, unsigned pos, char base, unsigned k)
{
    fhVal ^= rol(seedTab[(unsigned char)kmerSeq[pos]], k-1-pos);
    rhVal ^= rol(seedTab[(unsigned char)kmerSeq[pos]+cpOff], pos);
    kmerSeq[pos] = base;
    fhVal ^= rol(seedTab[(unsigned char)kmerSeq[pos]], k-1-pos);
    rhVal ^= rol(seedTab[(unsigned char)kmerSeq[pos]+cpOff], pos);
    return rhVal^fhVal;
}

/**
 * Compute multiple pseudo-independent hash values from a seed hash value.
 *
 * @param hashes array for storing computed hash values
 * @param seedVal seed value for multi-hash calculation
 * @param numHashes number of hash values to compute
 * @param k-mer size
 */
inline void multiHash(uint64_t hashes[], uint64_t seedVal, unsigned numHashes, unsigned k)
{
    for (unsigned i = 0; i < numHashes; i++) {
        hashes[i] = seedVal * (i ^ k * varSeed);
        hashes[i] ^= hashes[i] >> varShift;
    }
}

// spaced-seed hash values

/**
 * Calculate forward-strand spaced seed hash value of the base kmer, i.e. fhval(kmer_0)
 *
 * @param kVal set to forward-strand hash value for unmasked k-mer
 * @param seedSeq bitmask indicating "don't care" positions for hashing
 * @param kmerSeq k-mer to be hashed
 * @param k k-mer size
 * @return hash value for masked forward-strand k-mer
 */
inline uint64_t getFhval(uint64_t &kVal, const char * seedSeq, const char * kmerSeq, const unsigned k) {
    kVal=0;
    uint64_t sVal=0;
    for(unsigned i=0; i<k; i++) {
        kVal ^= rol(seedTab[(unsigned char)kmerSeq[i]], k-1-i);
        if(seedSeq[i]=='1')
            sVal ^= rol(seedTab[(unsigned char)kmerSeq[i]], k-1-i);
    }
    return sVal;
}

/**
 * Calculate reverse-strand spaced seed hash value of the base kmer, i.e. rhval(kmer_0)
 *
 * @param kVal set to reverse-strand hash value for unmasked k-mer
 * @param seedSeq bitmask indicating "don't care" positions for hashing
 * @param kmerSeq k-mer to be hashed
 * @param k k-mer size
 * @return hash for masked reverse-strand k-mer
 */
// reverse-strand spaced seed hash value of the base kmer, i.e. rhval(kmer_0)
inline uint64_t getRhval(uint64_t &kVal, const char * seedSeq, const char * kmerSeq, const unsigned k) {
    kVal=0;
    uint64_t sVal=0;
    for(unsigned i=0; i<k; i++) {
        kVal ^= rol(seedTab[(unsigned char)kmerSeq[i]+cpOff], i);
        if(seedSeq[i]=='1')
            sVal ^= rol(seedTab[(unsigned char)kmerSeq[i]+cpOff], i);
    }
    return sVal;
}

/**
 * Recursive forward-strand spaced seed hash value for next k-mer
 *
 * @param kVal hash value for current k-mer unmasked and in forward orientation
 * @param seedSeq bitmask indicating "don't care" positions for hashing
 * @param kmerSeq sequence for *current* k-mer (not the k-mer we are rolling into)
 * @param charIn new base we are rolling in from the right
 * @param k k-mer size
 * @return hash for masked k-mer in forward orientation
 */
inline uint64_t rollHashesRight(uint64_t &kVal, const char * seedSeq, const char * kmerSeq, const unsigned char charIn, const unsigned k) {
    const unsigned charOut = kmerSeq[0];
    kVal = rol(kVal, 1) ^ rol(seedTab[charOut], k) ^ seedTab[charIn];
    uint64_t sVal=kVal;
    for(unsigned i=1; i<k-1; i++) {
        if(seedSeq[i]!='1')
            sVal ^= rol(seedTab[(unsigned char)kmerSeq[i+1]], k-1-i);
    }
    return sVal;
}

/**
 * Recursive forward-strand spaced seed hash value for prev k-mer
 *
 * @param kVal hash value for current k-mer unmasked and in forward orientation
 * @param seedSeq bitmask indicating "don't care" positions for hashing
 * @param kmerSeq sequence for current k-mer (not the k-mer we are rolling into)
 * @param charIn new base we are rolling in from the left
 * @param k k-mer size
 * @return hash for masked k-mer in forward orientation
 */
inline uint64_t rollHashesLeft(uint64_t &kVal, const char * seedSeq, const char * kmerSeq, const unsigned char charIn, const unsigned k) {
    const unsigned charOut = kmerSeq[k-1];
    kVal = ror(kVal, 1) ^ ror(seedTab[charOut], 1) ^ rol(seedTab[charIn], k-1);
    uint64_t sVal=kVal;
    for(unsigned i=1; i<k-1; i++) {
        if(seedSeq[i]!='1')
            sVal ^= rol(seedTab[(unsigned char)kmerSeq[i-1]], k-1-i);
    }
    return sVal;
}

/**
 * Recursive canonical spaced seed hash value for next k-mer
 *
 * @param fkVal hash value for current k-mer unmasked and in forward orientation
 * @param rkVal hash value for current k-mer unmasked and in reverse complement orientation
 * @param seedSeq bitmask indicating "don't care" positions for hashing
 * @param kmerSeq sequence for current k-mer (not the k-mer we are rolling into)
 * @param charIn new base we are rolling in from the right
 * @param k k-mer size
 * @return canonical hash value for masked k-mer
 */
inline uint64_t rollHashesRight(uint64_t &fkVal, uint64_t &rkVal, const char * seedSeq, const char * kmerSeq, const unsigned char charIn, const unsigned k) {
    const unsigned charOut = kmerSeq[0];
    fkVal = rol(fkVal, 1) ^ rol(seedTab[charOut], k) ^ seedTab[charIn];
    rkVal = ror(rkVal, 1) ^ ror(seedTab[charOut+cpOff], 1) ^ rol(seedTab[charIn+cpOff], k-1);
    uint64_t fsVal=fkVal, rsVal=rkVal;
    for(unsigned i=1; i<k-1; i++) {
        if(seedSeq[i]!='1') {
            fsVal ^= rol(seedTab[(unsigned char)kmerSeq[i+1]], k-1-i);
            rsVal ^= rol(seedTab[(unsigned char)kmerSeq[i+1]+cpOff], i);
        }
    }
    return  rsVal^fsVal;
}

/**
 * Recursive canonical spaced seed hash value for prev k-mer
 *
 * @param fkVal hash value for current k-mer unmasked and in forward orientation
 * @param rkVal hash value for current k-mer unmasked and in reverse complement orientation
 * @param seedSeq bitmask indicating "don't care" positions for hashing
 * @param kmerSeq sequence for current k-mer (not the k-mer we are rolling into)
 * @param charIn new base we are rolling in from the left
 * @param k k-mer size
 * @return canonical hash value for masked k-mer
 */
inline uint64_t rollHashesLeft(uint64_t &fkVal, uint64_t &rkVal, const char * seedSeq, const char * kmerSeq, const unsigned char charIn, const unsigned k) {
    const unsigned charOut = kmerSeq[k-1];
    fkVal = ror(fkVal, 1) ^ ror(seedTab[charOut], 1) ^ rol(seedTab[charIn], k-1);
    rkVal = rol(rkVal, 1) ^ rol(seedTab[charOut+cpOff], k) ^ seedTab[charIn+cpOff];
    uint64_t fsVal=fkVal, rsVal=rkVal;
    for(unsigned i=1; i<k-1; i++) {
        if(seedSeq[i]!='1') {
            fsVal ^= rol(seedTab[(unsigned char)kmerSeq[i-1]], k-1-i);
            rsVal ^= rol(seedTab[(unsigned char)kmerSeq[i-1]+cpOff], i);
        }
    }
    return rsVal^fsVal;
}

/**
 * Change a single base and recompute spaced seed hash values
 *
 * @param fkVal hash value for current k-mer unmasked and in forward orientation
 * @param rkVal hash value for current k-mer unmasked and in reverse complement orientation
 * @param seedSeq bitmask indicating "don't care" positions for hashing
 * @param kmerSeq sequence for current k-mer
 * @param pos position of base to change
 * @param base new base value
 * @param k k-mer size
 * @return updated canonical hash value for masked k-mer
 */
inline uint64_t setBase(uint64_t& fkVal, uint64_t& rkVal, const char * seedSeq, char * kmerSeq, unsigned pos, char base, unsigned k)
{
    setBase(fkVal, rkVal, kmerSeq, pos, base, k);
    uint64_t fsVal=fkVal, rsVal=rkVal;
    for(unsigned i=0; i<k; i++) {
        if(seedSeq[i]!='1') {
            fsVal ^= rol(seedTab[(unsigned char)kmerSeq[i]], k-1-i);
            rsVal ^= rol(seedTab[(unsigned char)kmerSeq[i]+cpOff], i);
        }
    }
    return rsVal^fsVal;
}

#endif
