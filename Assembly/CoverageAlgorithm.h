#ifndef ASSEMBLY_COVERAGEALGORITHM_H
#define ASSEMBLY_COVERAGEALGORITHM_H 1

#include "Common/Histogram.h"
#include "Common/IOUtil.h"
#include "Common/Options.h" // for opt::rank
#include <fstream>

namespace AssemblyAlgorithms {

/** Return the k-mer coverage histogram. */
static inline
Histogram coverageHistogram(const SequenceCollectionHash& c)
{
	typedef SequenceCollectionHash Graph;

	Histogram h;
	for (Graph::const_iterator it = c.begin();
			it != c.end(); ++it) {
		if (it->second.deleted())
			continue;
		h.insert(it->second.getMultiplicity());
	}
	return h;
}

/** Calculate a k-mer coverage threshold from the given k-mer coverage
 * histogram. */
static inline
float calculateCoverageThreshold(const Histogram& h)
{
	float cov = h.firstLocalMinimum();
	if (opt::rank <= 0) {
		if (cov == 0)
			std::cout << "Unable to determine minimum k-mer coverage\n";
		else
			std::cout << "Minimum k-mer coverage is " << cov << std::endl;
	}

	for (unsigned iteration = 0; iteration < 100; iteration++) {
		Histogram trimmed = h.trimLow((unsigned)roundf(cov));
		if (opt::rank <= 0)
			logger(1) << "Coverage: " << cov << "\t"
				"Reconstruction: " << trimmed.size() << std::endl;

		unsigned median = trimmed.median();
		float cov1 = sqrt(median);
		if (cov1 == cov) {
			// The coverage threshold has converged.
			if (opt::rank <= 0)
				std::cout << "Using a coverage threshold of "
					<< (unsigned)roundf(cov) << "...\n"
					"The median k-mer coverage is " << median << "\n"
					"The reconstruction is " << trimmed.size()
					<< std::endl;
			if (!opt::db.empty()) {
				addToDb("coverageThreshold", (unsigned)roundf(cov));
				addToDb("medianKcoverage", median);
				addToDb("restruction", trimmed.size());
			}
			return cov;
		}
		cov = cov1;
	}
	if (opt::rank <= 0)
		std::cerr << "warning: coverage threshold did not converge"
			<< std::endl;
	return 0;
}

/** Set the coverage-related parameters e and c from the given k-mer
 * coverage histogram. */
static inline
void setCoverageParameters(const Histogram& h)
{
	if (!opt::coverageHistPath.empty() && opt::rank <= 0) {
		std::ofstream histFile(opt::coverageHistPath.c_str());
		assert_good(histFile, opt::coverageHistPath);
		histFile << h;
		assert(histFile.good());
	}

	float minCov = calculateCoverageThreshold(h);
	if (opt::rank <= 0) {
		if (minCov == 0)
			std::cout << "Unable to determine the "
				"k-mer coverage threshold" << std::endl;
		else
			std::cout << "The k-mer coverage threshold is " << minCov
				<< std::endl;
	}
	if (minCov < 2)
		minCov = 2;

	if ((int)opt::erode < 0) {
		opt::erode = (unsigned)roundf(minCov);
		if (opt::rank <= 0)
			std::cout << "Setting parameter e (erode) to "
				<< opt::erode << std::endl;
	}
	if ((int)opt::erodeStrand < 0) {
		opt::erodeStrand = minCov <= 2 ? 0 : 1;
		if (opt::rank <= 0)
			std::cout << "Setting parameter E (erodeStrand) to "
				<< opt::erodeStrand << std::endl;
	}
	if (opt::coverage < 0) {
		opt::coverage = minCov;
		if (opt::rank <= 0)
			std::cout << "Setting parameter c (coverage) to "
				<< opt::coverage << std::endl;
	}
}

/** Remove all k-mers with multiplicity lower than the given threshold */
static inline
size_t applyKmerCoverageThreshold(SequenceCollectionHash& c, unsigned kc)
{
	if (kc == 0)
		return 0;

	for (SequenceCollectionHash::iterator it = c.begin();
		it != c.end(); ++it) {
		if (it->second.getMultiplicity() < kc)
			it->second.setFlag(SF_DELETE);
	}

	return c.cleanup();
}

} // namespace AssemblyAlgorithms

#endif
