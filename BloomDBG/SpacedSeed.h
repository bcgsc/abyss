#ifndef SPACED_SEED_H
#define SPACED_SEED_H

#include <string>
#include <cassert>
#include <algorithm>

namespace SpacedSeed {

	/**
	 * Generate a spaced seed pattern (bitmask) for two equal-size
	 * k-mers separated by a gap.
	 *
	 * @param k width of spaced seed pattern
	 * @param K size of the individual k-mers. K must be <= k/2.
	 * @return spaced seed pattern for gapped k-mer pair
	 */
	static inline std::string kmerPair(unsigned k, unsigned K)
	{
		assert(K <= k/2);
		std::string seed(k, '0');
		std::fill(seed.begin(), seed.begin() + K, '1');
		std::fill(seed.rbegin(), seed.rbegin() + K, '1');
		return seed;
	}

	/**
	 * Generate a Quadratic Residue (QR) seed. The background theory
	 * for QR seeds is described in:
	 *
	 * Egidi, Lavinia, and Giovanni Manzini. "Multiple seeds
	 * sensitivity using a single seed with threshold." Journal of
	 * bioinformatics and computational biology 13.04 (2015): 1550011.
	 *
	 * @param len desired length of QR seed. `len` must
	 * be prime and >= 11.
	 * @return a QR seed represented as a std::string
	 * of 0's and 1's
	 */
	static inline std::string qrSeed(unsigned len)
	{
		assert(len >= 11);
		std::string seed(len, '1');
		for (size_t i = 0; i < len; ++i) {
			for (size_t j = 1; j < len; ++j) {
				if (j*j % len == i) {
					seed.at(i) = '0';
					break;
				}
			}
		}
		return seed;
	}

	/**
	 * Generate a spaced seed pattern (bitmask) for two equal-length
	 * Quadratic Residue (QR) seeds separated by a gap.  The first
	 * QR seed is in the usual orientation and the second QR is reversed,
	 * so that the overall pattern is symmetric.
	 *
	 * @param k width of the spaced seed pattern
	 * @param qrSeedLen width of the individual QR seeds.
	 * qrSeedLen must be a prime number >= 11 and must also be <= k/2.
	 * @return spaced seed pattern for gapped QR seed pair
	 */
	static inline std::string qrSeedPair(unsigned k, unsigned qrSeedLen)
	{
		assert(qrSeedLen <= k/2);
		std::string seed(k, '0');
		std::string qrSeed = SpacedSeed::qrSeed(qrSeedLen);
		std::copy(qrSeed.begin(), qrSeed.end(), seed.begin());
		std::reverse(qrSeed.begin(), qrSeed.end());
		std::copy(qrSeed.rbegin(), qrSeed.rend(), seed.rbegin());
		return seed;
	}

}

#endif
