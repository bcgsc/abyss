#ifndef RRESOLVER_RALGORITHMS_SHORT_H
#define RRESOLVER_RALGORITHMS_SHORT_H 1

#include "BloomFilters.h"
#include "Common/Histogram.h"
#include "Contigs.h"

#include <cassert>
#include <climits>
#include <limits>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>

const int MIN_MARGIN = 2;
const int R_HEURISTIC = 45;
const double R_HEURISTIC_A = 0.49;
const double R_HEURISTIC_B = 63.5;
const int MAX_SUBITERATIONS = 2;
const long HIST_SAMPLE_SIZE = LONG_MAX;
const long REPEAT_CASES_LIMIT = LONG_MAX;
const long READ_STATS_SAMPLE_SIZE = 100000;
const double READ_BATCH_FRACTION_THRESHOLD = 0.15;
const long PATH_COMBINATIONS_MULTITHREAD_THRESHOLD = 5000;
const double SUPPORTED_PATHS_MIN = 0.15;
const double COV_APPROX_FORMULA_FACTOR = 4.00;
const double SPACED_SEEDS_SNP_FRACTION = 1.00;
const int MAX_READ_SIZE = 300;

namespace opt {

/** Read Bloom filter size in bytes. */
extern size_t bloomSize;

/** The number of parallel threads */
extern int threads;

/** Prefix for the histogram files */
extern std::string histPrefix;

/** Name of the file to write resulting graph to */
extern std::string outputGraphPath;

/** Name of the file to write resulting graph to */
extern std::string outputContigsPath;

/** Number of kmers required to be found for a path to be supported */
extern int threshold;

/** Number of Rmers to extract per read */
extern int extract;

/** Minimum number of sliding window moves */
extern int minTests;

/** Maximum number of sliding window moves */
extern int maxTests;

/** Maximum number of branching paths */
extern int branching;

/** Explicitly specified r values. */
extern std::vector<int> rValues;

/** Explicitly set coverage approximation factors. */
extern std::vector<double> covApproxFactors;

/** Base pair quality threshold that large k-mers bases should have at minimum. */
extern int readQualityThreshold;

/** Flag indicating whether error correction is enabled */
extern int errorCorrection;

/** Name of the file to write supported paths to */
extern std::string outputSupportedPathsPath;

/** Name of the file to write unsupported paths to */
extern std::string outputUnsupportedPathsPath;

extern unsigned k; // used by ContigProperties

extern int format; // used by ContigProperties

}

class Support
{
  public:

	enum class UnknownReason : uint8_t {
		UNDETERMINED = 0, // Not yet processed
		TOO_MANY_COMBINATIONS, // Branching out exploded beyond a threshold
		OVER_MAX_TESTS, // Planned tests was above the threshold
		POSSIBLE_TESTS_LT_PLANNED, // Planned tests could not be carried out due to low coverage
		WINDOW_NOT_LONG_ENOUGH, // Window too small / repeat too large to all the planned tests
		HEAD_SHORTER_THAN_MARGIN, // One of the branches to the left was too short for planned tests
		TAIL_SHORTER_THAN_MARGIN, // One of the branches to the right was too short for planned tests
		DIFFERENT_CULPRIT // The path was fine, but another path in this repeat was unknown so all the paths
		                  // for this repeat became unknown
	};

	Support() = default;

	Support(UnknownReason unknownReason): unknownReason(unknownReason) {};

	Support(int calculatedTests, UnknownReason unknownReason)
	  : found(-1)
	  , tests(-1)
		, unknownReason(unknownReason)
	{
		assert(calculatedTests >= 0);
		if (calculatedTests > std::numeric_limits<decltype(this->calculatedTests)>::max()) { this->calculatedTests = std::numeric_limits<decltype(this->calculatedTests)>::max(); } else {
			this->calculatedTests = calculatedTests;
		}
	}

	Support(int8_t found, int8_t tests)
	  : found(found)
	  , tests(tests)
	{
		assert(found >= 0);
		assert(tests > 0);
		assert(calculatedTests == -1);
	}

	Support(int8_t found, int8_t tests, int calculatedTests)
	  : found(found)
	  , tests(tests)
	{
		assert(found >= 0);
		assert(tests > 0);
		assert(calculatedTests >= 0);
		if (calculatedTests > std::numeric_limits<decltype(this->calculatedTests)>::max()) { this->calculatedTests = std::numeric_limits<decltype(this->calculatedTests)>::max(); } else {
			this->calculatedTests = calculatedTests;
		}
	}

	bool unknown() const { return tests == -1; }

	void reset()
	{
		found = -1;
		tests = -1;
	}

	bool good() const { return unknown() || found >= opt::threshold; }

	bool operator>(const Support& s) { return found > s.found; }

	bool operator<(const Support& s) { return found < s.found; }

	int8_t found = -1;
	int8_t tests = -1;
	int8_t calculatedTests = -1;
	UnknownReason unknownReason = UnknownReason::UNDETERMINED;
};

class ReadSize
{

  public:
	ReadSize(int size)
	  : size(size)
	{}

	double getFractionOfTotal() const { return double(sampleCount) / double(readsSampleSize); }

	int size;
	std::vector<int> rValues;
	Histogram qualThresholdPositions;
	long sampleCount = 0;
	double covApproxFactor = COV_APPROX_FORMULA_FACTOR;

	static long readsSampleSize;
	static std::vector<ReadSize> readSizes;
	static ReadSize current;
};

void
resolveShort(
    const std::vector<std::string>& readFilepaths,
    ImaginaryContigPaths& supportedPaths,
    ImaginaryContigPaths& unsupportedPaths);

#endif
