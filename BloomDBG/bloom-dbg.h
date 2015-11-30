#ifndef BLOOM_DBG_H
#define BLOOM_DBG_H 1

#include <cmath>

/**
 * Parameters for creating a Bloom filter.
 */
struct BloomParams
{
	/** Bloom filter size in bits */
	std::size_t size;
	/** number of hash functions */
	unsigned hashes;
};

/**
 * Calculate optimal Bloom filter params based on
 * genome size, read coverage, sequencing error
 * rate, k-mer size, and target false positive rate (FPR).
 *
 * @param genomeSize genome size in bp
 * @param readCoverage fold read coverage
 * @param seqErrorRate sequencing error rate between 0 and 1
 * @param k k-mer length in bp
 * @param fpr target false positive rate
 *
 * @return BloomParams struct containing Bloom filter size
 * and number of hash funcitons
 */
BloomParams calcBloomParams(std::size_t genomeSize, double readCoverage,
	double seqErrorRate, unsigned k, double fpr)
{
	BloomParams params;

	/* approx number of k-mers that will be inserted into Bloom filter */
	std::size_t n = genomeSize + genomeSize * readCoverage * seqErrorRate * k;

	/* see https://en.wikipedia.org/wiki/Bloom_filter */
	params.size = (std::size_t)(-(double)n*log(fpr)/pow(log(2.0), 2));
	params.hashes = (unsigned)ceil((double)params.size*log(2)/n);

	return params;
}

#endif
