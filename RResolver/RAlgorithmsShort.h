#ifndef RRESOLVER_RALGORITHMS_SHORT_H
#define RRESOLVER_RALGORITHMS_SHORT_H 1

#include "RUtils.h"
#include "BloomFilters.h"
#include "Contigs.h"
#include "SequenceTree.h"
#include "Common/Histogram.h"

#include <vector>
#include <string>
#include <cassert>
#include <fstream>
#include <climits>
#include <cmath>

const int    MIN_MARGIN = 2;
const int    MAX_TESTS_OFFSET = 16;
const int    R_VALUES_STEP = 20;
const int    R_STEPS_MAX = 1;
const int    R_MAX_K_DIFF = 80;
const int    MAX_SUBITERATIONS = 2;
const long   HIST_SAMPLE_SIZE = LONG_MAX;
const long   REPEAT_CASES_LIMIT = LONG_MAX;
const int    RMER_QUALITY_THRESHOLD = 30;
const long   READ_STATS_SAMPLE_SIZE = 100000;
const double READ_BATCH_FRACTION_THRESHOLD = 0.30;
const long   PATH_COMBINATIONS_MULTITHREAD_THRESHOLD = 5000;
const double SUPPORTED_PATHS_MIN = 0.15;
const double COV_APPROX_FORMULA_FACTOR = 5.00;
const double SPACED_SEEDS_SNP_FRACTION = 1.00;
const int    MAX_READ_SIZE = 300;

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

  /** Minimum number of sliding window moves */
  extern int min_tests;

  /** Maximum number of branching paths */
  extern int branching;

  /** Flag indicating whether error correction is enabled */
  extern int error_correction;

  /** Name of the file to write supported paths to */
  extern std::string outputSupportedPathsPath;

  /** Name of the file to write unsupported paths to */
  extern std::string outputUnsupportedPathsPath;

  extern unsigned k;  // used by ContigProperties

  extern int format;  // used by ContigProperties

}

class Support {

public:

  Support() = default;

  Support(int calculatedTests): found(-1), tests(-1), calculatedTests(calculatedTests) {
    assert(calculatedTests >= 0);
  }

  Support(int found, int tests):
    found(found), tests(tests)
  {
    assert(found >= 0);
    assert(tests > 0);
    assert(calculatedTests == -1);
  }

  Support(int found, int tests, int calculatedTests):
    found(found), tests(tests), calculatedTests(calculatedTests)
  {
    assert(found >= 0);
    assert(tests > 0);
    assert(calculatedTests >= 0);
  }

  bool unknown() const {
    return tests == -1;
  }

  void reset() {
    found = -1;
    tests = -1;
  }

  bool good() const {
    return unknown() || found >= opt::threshold;
  }

  bool operator>(const Support& s) {
    return found > s.found;
  }

  bool operator<(const Support& s) {
    return found < s.found;
  }

  int found = -1;
  int tests = -1;
  int calculatedTests = -1;

};

class ReadBatch {

public:

  ReadBatch(int size): size(size) {}

  double getFractionOfTotal() const { return double(sampleCount) / double(readsSampleSize); }

  int size;
  std::vector<int> rValues;
  Histogram qualThresholdPositions;
  long sampleCount = 0;

  static long readsSampleSize;
  static std::vector<ReadBatch> batches;
  static ReadBatch current;

};

void
resolveShort(const std::vector<std::string>& readFilepaths,
             ImaginaryContigPaths& supportedPaths,
             ImaginaryContigPaths& unsupportedPaths);

#endif
