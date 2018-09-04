/*
 *
 * nthash.hpp
 * Author: Hamid Mohamadi
 * Genome Sciences Centre,
 * British Columbia Cancer Agency
 */

#ifndef NT_HASH_H
#define NT_HASH_H

#include <stdint.h>

// offset for the complement base in the random seeds table
const uint8_t cpOff = 0x07;

// shift for gerenerating multiple hash values
const int multiShift = 27;

// seed for gerenerating multiple hash values
static const uint64_t multiSeed = 0x90b45d39fb6da1fa;

// 64-bit random seeds corresponding to bases and their complements
static const uint64_t seedA = 0x3c8bfbb395c60474;
static const uint64_t seedC = 0x3193c18562a02b4c;
static const uint64_t seedG = 0x20323ed082572324;
static const uint64_t seedT = 0x295549f54be24456;
static const uint64_t seedN = 0x0000000000000000;

static const uint64_t seedTab[256] = {
    seedN, seedT, seedN, seedG, seedA, seedN, seedN, seedC, // 0..7
    seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 8..15
    seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 16..23
    seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 24..31
    seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 32..39
    seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 40..47
    seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 48..55
    seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 56..63
    seedN, seedA, seedN, seedC, seedN, seedN, seedN, seedG, // 64..71
    seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 72..79
    seedN, seedN, seedN, seedN, seedT, seedN, seedN, seedN, // 80..87
    seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 88..95
    seedN, seedA, seedN, seedC, seedN, seedN, seedN, seedG, // 96..103
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

// rotate "v" to the left 1 position
inline uint64_t rol1(const uint64_t v) {
    return (v << 1) | (v >> 63);
}

// rotate "v" to the right by 1 position
inline uint64_t ror1(const uint64_t v) {
    return (v >> 1) | (v << 63);
}

// rotate 31-left bits of "v" to the left by "s" positions
inline uint64_t rol31(const uint64_t v, unsigned s) {
    s%=31;
    return ((v << s) | (v >> (31 - s))) & 0x7FFFFFFF;
}

// rotate 33-right bits of "v" to the left by "s" positions
inline uint64_t rol33(const uint64_t v, unsigned s) {
    s%=33;
    return ((v << s) | (v >> (33 - s))) & 0x1FFFFFFFF;
}

// swap bit 0 with bit 33 in "v"
inline uint64_t swapbits033(const uint64_t v) {
    uint64_t x = (v ^ (v >> 33)) & 1;
    return v ^ (x | (x << 33));
}

// swap bit 32 with bit 63 in "v"
inline uint64_t swapbits3263(const uint64_t v) {
    uint64_t x = ((v >> 32) ^ (v >> 63)) & 1;
    return v ^ ((x << 32) | (x << 63));
}

// forward-strand hash value of the base kmer, i.e. fhval(kmer_0)
inline uint64_t NTF64(const char * kmerSeq, const unsigned k) {
    uint64_t hVal=0;
    for(unsigned i=0; i<k; i++) {
        hVal = rol1(hVal);
        hVal = swapbits033(hVal);
        hVal ^= seedTab[(unsigned char)kmerSeq[i]];
    }
    return hVal;
}

// reverse-strand hash value of the base kmer, i.e. rhval(kmer_0)
inline uint64_t NTR64(const char * kmerSeq, const unsigned k) {
    uint64_t hVal=0;
    for(unsigned i=0; i<k; i++) {
        hVal = rol1(hVal);
        hVal = swapbits033(hVal);
        hVal ^= seedTab[(unsigned char)kmerSeq[k-1-i]&cpOff];
    }
    return hVal;
}

// forward-strand ntHash for sliding k-mers
inline uint64_t NTF64(const uint64_t fhVal, const unsigned k, const unsigned char charOut, const unsigned char charIn) {
    uint64_t hVal = rol1(fhVal);
    hVal = swapbits033(hVal);
    hVal ^= seedTab[charIn];
    uint64_t lBits = seedTab[charOut] >> 33;
    uint64_t rBits = seedTab[charOut] & 0x1FFFFFFFF;
    uint64_t sOut = (rol31(lBits,k) << 33) | (rol33(rBits,k));
    hVal ^= sOut;
    return hVal;
}

// reverse-complement ntHash for sliding k-mers
inline uint64_t NTR64(const uint64_t rhVal, const unsigned k, const unsigned char charOut, const unsigned char charIn) {
    uint64_t lBits = seedTab[charIn&cpOff] >> 33;
    uint64_t rBits = seedTab[charIn&cpOff] & 0x1FFFFFFFF;
    uint64_t sIn = (rol31(lBits,k) << 33) | (rol33(rBits,k));
    uint64_t hVal = rhVal ^ sIn;
    hVal ^= seedTab[charOut&cpOff];
    hVal = ror1(hVal);
    hVal = swapbits3263(hVal);
    return hVal;
}

// canonical ntBase
inline uint64_t NTC64(const char * kmerSeq, const unsigned k) {
    uint64_t fhVal=0, rhVal=0;
    fhVal=NTF64(kmerSeq, k);
    rhVal=NTR64(kmerSeq, k);
    return (rhVal<fhVal)? rhVal : fhVal;
}

// canonical ntHash
inline uint64_t NTC64(const char * kmerSeq, const unsigned k, uint64_t& fhVal, uint64_t& rhVal) {
    fhVal = NTF64(kmerSeq, k);
    rhVal = NTR64(kmerSeq, k);
    return (rhVal<fhVal)? rhVal : fhVal;
}

// canonical ntHash for sliding k-mers
inline uint64_t NTC64(const unsigned char charOut, const unsigned char charIn, const unsigned k, uint64_t& fhVal, uint64_t& rhVal) {
    fhVal = NTF64(fhVal, k, charOut, charIn);
    rhVal = NTR64(rhVal, k, charOut, charIn);
    return (rhVal<fhVal)? rhVal : fhVal;
}

// forward-strand ntHash for sliding k-mers to the left
inline uint64_t NTF64L(const uint64_t rhVal, const unsigned k, const unsigned char charOut, const unsigned char charIn) {
    uint64_t lBits = seedTab[charIn] >> 33;
    uint64_t rBits = seedTab[charIn] & 0x1FFFFFFFF;
    uint64_t sIn = (rol31(lBits,k) << 33) | (rol33(rBits,k));
    uint64_t hVal = rhVal ^ sIn;
    hVal ^= seedTab[charOut];
    hVal = ror1(hVal);
    hVal = swapbits3263(hVal);
    return hVal;
}

// reverse-complement ntHash for sliding k-mers to the left
inline uint64_t NTR64L(const uint64_t fhVal, const unsigned k, const unsigned char charOut, const unsigned char charIn) {
    uint64_t hVal = rol1(fhVal);
    hVal = swapbits033(hVal);
    hVal ^= seedTab[charIn&cpOff];
    uint64_t lBits = seedTab[charOut&cpOff] >> 33;
    uint64_t rBits = seedTab[charOut&cpOff] & 0x1FFFFFFFF;
    uint64_t sOut = (rol31(lBits,k) << 33) | (rol33(rBits,k));
    hVal ^= sOut;
    return hVal;
}

// canonical ntHash for sliding k-mers to the left
inline uint64_t NTC64L(const unsigned char charOut, const unsigned char charIn, const unsigned k, uint64_t& fhVal, uint64_t& rhVal) {
    fhVal = NTF64L(fhVal, k, charOut, charIn);
    rhVal = NTR64L(rhVal, k, charOut, charIn);
    return (rhVal<fhVal)? rhVal : fhVal;
}

// ntBase with seeding option
inline uint64_t NTF64(const char * kmerSeq, const unsigned k, const unsigned seed) {
    uint64_t hVal=NTF64(kmerSeq, k);
    if(seed==0) return hVal;
    hVal *= seed ^ k * multiSeed;
    hVal ^= hVal >> multiShift;
    return hVal;
}

// canonical ntBase with seeding option
inline uint64_t NTC64(const char * kmerSeq, const unsigned k, const unsigned seed) {
    uint64_t hVal = NTC64(kmerSeq,k);
    if(seed==0) return hVal;
    hVal *= seed ^ k * multiSeed;
    hVal ^= hVal >> multiShift;
    return hVal;
}

// multihash ntHash, ntBase
inline void NTM64(const char * kmerSeq, const unsigned k, const unsigned m, uint64_t *hVal) {
    uint64_t bVal=0, tVal=0;
    bVal = NTF64(kmerSeq, k);
    hVal[0] = bVal;
    for(unsigned i=1; i<m; i++) {
        tVal = bVal * (i ^ k * multiSeed);
        tVal ^= tVal >> multiShift;
        hVal[i] =  tVal;
    }
}

// one extra hash for given base hash
inline uint64_t NTE64(const uint64_t hVal, const unsigned k, const unsigned i) {
    uint64_t tVal = hVal;
    tVal *= (i ^ k * multiSeed);
    tVal ^= tVal >> multiShift;
    return tVal;
}

// multihash ntHash for sliding k-mers
inline void NTM64(const unsigned char charOut, const unsigned char charIn, const unsigned k, const unsigned m, uint64_t *hVal) {
    uint64_t bVal=0, tVal=0;
    bVal = NTF64(hVal[0], k, charOut, charIn);
    hVal[0] = bVal;
    for(unsigned i=1; i<m; i++) {
        tVal = bVal * (i ^ k * multiSeed);
        tVal ^= tVal >> multiShift;
        hVal[i] =  tVal;
    }
}

// canonical multihash ntBase
inline void NTMC64(const char * kmerSeq, const unsigned k, const unsigned m, uint64_t *hVal) {
    uint64_t bVal=0, tVal=0;
    bVal = NTC64(kmerSeq, k);
    hVal[0] = bVal;
    for(unsigned i=1; i<m; i++) {
        tVal = bVal * (i ^ k * multiSeed);
        tVal ^= tVal >> multiShift;
        hVal[i] =  tVal;
    }
}

// canonical multihash ntHash
inline void NTMC64(const char * kmerSeq, const unsigned k, const unsigned m, uint64_t& fhVal, uint64_t& rhVal, uint64_t *hVal) {
    uint64_t bVal=0, tVal=0;
    bVal = NTC64(kmerSeq, k, fhVal, rhVal);
    hVal[0] = bVal;
    for(unsigned i=1; i<m; i++) {
        tVal = bVal * (i ^ k * multiSeed);
        tVal ^= tVal >> multiShift;
        hVal[i] =  tVal;
    }
}

// canonical multihash ntHash for sliding k-mers
inline void NTMC64(const unsigned char charOut, const unsigned char charIn, const unsigned k, const unsigned m, uint64_t& fhVal, uint64_t& rhVal, uint64_t *hVal) {
    uint64_t bVal=0, tVal=0;
    bVal = NTC64(charOut, charIn, k, fhVal, rhVal);
    hVal[0] = bVal;
    for(unsigned i=1; i<m; i++) {
        tVal = bVal * (i ^ k * multiSeed);
        tVal ^= tVal >> multiShift;
        hVal[i] =  tVal;
    }
}

/*
 * ignoring k-mers containing nonACGT using ntHash function
*/

// canonical ntBase
inline bool NTC4(const char *kmerSeq, const unsigned k, uint64_t& hVal, unsigned& locN) {
    hVal=0;
    locN=0;
    uint64_t fhVal=0,rhVal=0;
    for(int i=k-1; i>=0; i--) {
        if(seedTab[(unsigned char)kmerSeq[i]]==seedN) {
            locN=i;
            return false;
        }
        fhVal = rol1(fhVal);
        fhVal = swapbits033(fhVal);
        fhVal ^= seedTab[(unsigned char)kmerSeq[k-1-i]];

        rhVal = rol1(rhVal);
        rhVal = swapbits033(rhVal);
        rhVal ^= seedTab[(unsigned char)kmerSeq[i]&cpOff];
    }
    hVal = (rhVal<fhVal)? rhVal : fhVal;
    return true;
}

// canonical multihash ntBase
inline bool NTMC64(const char *kmerSeq, const unsigned k, const unsigned m, unsigned& locN, uint64_t* hVal) {
    uint64_t bVal=0, tVal=0, fhVal=0, rhVal=0;
    locN=0;
    for(int i=k-1; i>=0; i--) {
        if(seedTab[(unsigned char)kmerSeq[i]]==seedN) {
            locN=i;
            return false;
        }
        fhVal = rol1(fhVal);
        fhVal = swapbits033(fhVal);
        fhVal ^= seedTab[(unsigned char)kmerSeq[k-1-i]];

        rhVal = rol1(rhVal);
        rhVal = swapbits033(rhVal);
        rhVal ^= seedTab[(unsigned char)kmerSeq[i]&cpOff];
    }
    bVal = (rhVal<fhVal)? rhVal : fhVal;
    hVal[0] = bVal;
    for(unsigned i=1; i<m; i++) {
        tVal = bVal * (i ^ k * multiSeed);
        tVal ^= tVal >> multiShift;
        hVal[i] =  tVal;
    }
    return true;
}

// canonical ntHash
inline bool NTC64(const char *kmerSeq, const unsigned k, uint64_t& fhVal, uint64_t& rhVal, uint64_t& hVal, unsigned& locN) {
    hVal=fhVal=rhVal=0;
    locN=0;
    for(int i=k-1; i>=0; i--) {
        if(seedTab[(unsigned char)kmerSeq[i]]==seedN) {
            locN=i;
            return false;
        }
        fhVal = rol1(fhVal);
        fhVal = swapbits033(fhVal);
        fhVal ^= seedTab[(unsigned char)kmerSeq[k-1-i]];

        rhVal = rol1(rhVal);
        rhVal = swapbits033(rhVal);
        rhVal ^= seedTab[(unsigned char)kmerSeq[i]&cpOff];
    }
    hVal = (rhVal<fhVal)? rhVal : fhVal;
    return true;
}

// canonical multihash ntHash
inline bool NTMC64(const char *kmerSeq, const unsigned k, const unsigned m, uint64_t& fhVal, uint64_t& rhVal, unsigned& locN, uint64_t* hVal) {
    fhVal=rhVal=0;
    uint64_t bVal=0, tVal=0;
    locN=0;
    for(int i=k-1; i>=0; i--) {
        if(seedTab[(unsigned char)kmerSeq[i]]==seedN) {
            locN=i;
            return false;
        }
        fhVal = rol1(fhVal);
        fhVal = swapbits033(fhVal);
        fhVal ^= seedTab[(unsigned char)kmerSeq[k-1-i]];

        rhVal = rol1(rhVal);
        rhVal = swapbits033(rhVal);
        rhVal ^= seedTab[(unsigned char)kmerSeq[i]&cpOff];
    }
    bVal = (rhVal<fhVal)? rhVal : fhVal;
    hVal[0] = bVal;
    for(unsigned i=1; i<m; i++) {
        tVal = bVal * (i ^ k * multiSeed);
        tVal ^= tVal >> multiShift;
        hVal[i] =  tVal;
    }
    return true;
}

// strand-aware canonical multihash ntHash
inline bool NTMC64(const char *kmerSeq, const unsigned k, const unsigned m, uint64_t& fhVal, uint64_t& rhVal, unsigned& locN, uint64_t* hVal, bool& hStn) {
    fhVal=rhVal=0;
    uint64_t bVal=0, tVal=0;
    locN=0;
    for(int i=k-1; i>=0; i--) {
        if(seedTab[(unsigned char)kmerSeq[i]]==seedN) {
            locN=i;
            return false;
        }
        fhVal = rol1(fhVal);
        fhVal = swapbits033(fhVal);
        fhVal ^= seedTab[(unsigned char)kmerSeq[k-1-i]];
        
        rhVal = rol1(rhVal);
        rhVal = swapbits033(rhVal);
        rhVal ^= seedTab[(unsigned char)kmerSeq[i]&cpOff];
    }
    hStn = rhVal<fhVal;
    bVal = hStn? rhVal : fhVal;
    hVal[0] = bVal;
    for(unsigned i=1; i<m; i++) {
        tVal = bVal * (i ^ k * multiSeed);
        tVal ^= tVal >> multiShift;
        hVal[i] =  tVal;
    }
    return true;
}

// starnd-aware canonical multihash ntHash for sliding k-mers
inline void NTMC64(const unsigned char charOut, const unsigned char charIn, const unsigned k, const unsigned m, uint64_t& fhVal, uint64_t& rhVal, uint64_t *hVal, bool &hStn) {
    uint64_t bVal=0, tVal=0;
    bVal = NTC64(charOut, charIn, k, fhVal, rhVal);
    hStn = rhVal<fhVal;
    hVal[0] = bVal;
    for(unsigned i=1; i<m; i++) {
        tVal = bVal * (i ^ k * multiSeed);
        tVal ^= tVal >> multiShift;
        hVal[i] =  tVal;
    }
}

// masking canonical ntHash using spaced seed pattern
inline uint64_t maskHash(uint64_t &fkVal, uint64_t &rkVal, const char * seedSeq, const char * kmerSeq, const unsigned k) {
    uint64_t fsVal=fkVal, rsVal=rkVal;
    for(unsigned i=0; i<k; i++) {
        if(seedSeq[i]!='1') {
            uint64_t lfBits = seedTab[(unsigned char)kmerSeq[i]] >> 33;
            uint64_t rfBits = seedTab[(unsigned char)kmerSeq[i]] & 0x1FFFFFFFF;
            uint64_t sfMask = (rol31(lfBits,k-1-i) << 33) | (rol33(rfBits,k-1-i));
            fsVal ^= sfMask;

            uint64_t lrBits = seedTab[(unsigned char)kmerSeq[i]&cpOff] >> 33;
            uint64_t rrBits = seedTab[(unsigned char)kmerSeq[i]&cpOff] & 0x1FFFFFFFF;
            uint64_t srMask = (rol31(lrBits,i) << 33) | (rol33(rrBits,i));
            rsVal ^= srMask;
        }
    }
    return (rsVal<fsVal)? rsVal : fsVal;
}

#endif
