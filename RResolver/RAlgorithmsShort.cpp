#include "RAlgorithmsShort.h"
#include "vendor/nthash/stHashIterator.hpp"
#include "vendor/nthash/ntHashIterator.hpp"

#include "btllib/include/btllib/seq_reader.hpp"

#include <omp.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include <algorithm>

static unsigned char BASES[] = { 'A', 'C', 'T', 'G' };

typedef std::map<long, std::map<long, Support>> SupportMap;
typedef std::map<long, SupportMap> RepeatSupportMap;

long ReadBatch::readsSampleSize = 0;
std::vector<ReadBatch> ReadBatch::batches;
ReadBatch ReadBatch::current(0);

class FractionHistogram : public Histogram {

public:

  void insert(double fraction) {
    assert(fraction >= 0);
    assert(fraction <= 1);
    Histogram::insert(int(fraction * 100));
  }

  friend std::ostream& operator<<(std::ostream& o, const FractionHistogram& h)
	{
		o << (Histogram&)h;
    if (h.size() == 0 || (--h.end())->first != 100) {
      o << 100 << "\t0\n";
    }
		return o;
	}

};

class Resolution {

public:

  Resolution(const ReadBatch& batch, int r): batch(batch), r(r) {}

  RepeatSupportMap repeatSupportMap;
  const ReadBatch& batch;
  int r;
  Histogram findsHistogram;
  FractionHistogram fractionFindsHistogram;
  Histogram calculatedTestsHistogram;
  bool failed = false;

};

static int getMinWindowLength(const int tests, const int repeatSize, const int minMargin) {
  return tests - 1 + minMargin + repeatSize + minMargin;
}

static bool windowLongEnough(const int windowSize, const int tests, const int repeatSize, const int minMargin) {
  return windowSize >= getMinWindowLength(tests, repeatSize, minMargin);
}

//static int numOfTests(const int repeatSize, const int windowSize, const int minMargin) {
//  return windowSize - minMargin - repeatSize - minMargin + 1;
//}

static int getMargin(const int windowSize, const int tests, const int repeatSize, const int minMargin) {
  assert(windowLongEnough(windowSize, tests, repeatSize, minMargin));
  const int requiredSeqSize = windowSize + tests - 1;
  const int margin = (requiredSeqSize - repeatSize + 1) / 2;
  assert(margin >= minMargin);
  return margin;
}

static bool determineShortReadStats(const std::vector<std::string>& readFilenames) {
  if (opt::verbose) {
    std::cerr << "Determining read stats..." << std::endl;
  }
  ReadBatch::batches.clear();
  #pragma omp parallel
  #pragma omp single
  {
    for (const auto filename : readFilenames)
    #pragma omp task firstprivate(filename)
    {
      Histogram hist;
      std::map<int, Histogram> qualThresholdPositionsHists;

      btllib::SeqReader reader(filename);
      for (btllib::SeqReader::Record record; (record = reader.read()) && (record.num < READ_STATS_SAMPLE_SIZE);) {
        if (record.seq.size() > MAX_READ_SIZE) { continue; }
        hist.insert(record.seq.size());
        for (int j = record.qual.size() - 1; j >= 0; j--) {
          if (record.qual[j] >= RMER_QUALITY_THRESHOLD) {
            qualThresholdPositionsHists[record.seq.size()].insert(j);
            break;
          }
        }
      }

      #pragma omp critical (ReadBatches)
      {
        for (const auto& i : hist) {
          ReadBatch* batch = nullptr;
          bool found = false;
          for (auto& b : ReadBatch::batches) {
            if (b.size == i.first) {
              found = true;
              batch = &b;
              break;
            }
          }
          if (!found) {
            ReadBatch::batches.push_back(ReadBatch(i.first));
            batch = &(ReadBatch::batches.back());
          }
          batch->sampleCount += i.second;
          auto& qualHist = batch->qualThresholdPositions;
          for (const auto& q: qualThresholdPositionsHists[i.first]) {
            for (size_t n = 0; n < q.second; n++) {
              qualHist.insert(q.first);
            }
          }
        }
      }
    }
    #pragma omp taskwait
  }

  ReadBatch::readsSampleSize = 0;
  for (const auto& batch : ReadBatch::batches) {
    ReadBatch::readsSampleSize += batch.sampleCount;
  }

  std::sort(ReadBatch::batches.begin(), ReadBatch::batches.end(), [] (ReadBatch a, ReadBatch b) { return a.sampleCount > b.sampleCount; });

  if (ReadBatch::batches.size() == 0) {
    std::cerr << "Insufficient number of short reads. Finishing..." << std::endl;
    return false;
  }

  if (ReadBatch::batches[0].getFractionOfTotal() < READ_BATCH_FRACTION_THRESHOLD) {
    std::cerr << "Insufficient reads of same size. Finishing..." << std::endl;
    return false;
  }

  std::vector<ReadBatch> batchesFiltered;
  for (const auto& b: ReadBatch::batches) {
    if (b.getFractionOfTotal() >= READ_BATCH_FRACTION_THRESHOLD) {
      batchesFiltered.push_back(b);
    }
  }
  ReadBatch::batches = batchesFiltered;

  std::sort(ReadBatch::batches.begin(), ReadBatch::batches.end(), [] (ReadBatch a, ReadBatch b) { return a.size < b.size; });
  if (opt::verbose) {
    std::cerr << "Read lengths determined to be: " << std::fixed;
    std::cerr << ReadBatch::batches[0].size << " (" << (ReadBatch::batches[0].getFractionOfTotal() * 100.0) << "%)";
    for (size_t i = 1; i < ReadBatch::batches.size(); i++) {
      std::cerr << ", " << ReadBatch::batches[i].size << " (" << (ReadBatch::batches[i].getFractionOfTotal() * 100.0) << "%)";
    }
    std::cerr << std::defaultfloat << std::endl;
  }
  
  for (size_t i = 0; i < ReadBatch::batches.size(); i++) {
    auto& batch = ReadBatch::batches[i];
    int r = batch.qualThresholdPositions.median() - opt::threshold + 1;
    int prevR = opt::k;
    if (i > 0) {
      prevR = ReadBatch::batches[i - 1].rValues.back();
    }
    int steps = 0;
    while ((r - prevR > R_VALUES_STEP) && (steps < R_STEPS_MAX)) {
      if (r - opt::k > R_MAX_K_DIFF) { r = opt::k + R_MAX_K_DIFF; }
      batch.rValues.push_back(r);
      r -= R_VALUES_STEP;
      steps++;
    }
    std::reverse(batch.rValues.begin(), batch.rValues.end());
  }

  if (opt::verbose) {
    std::cerr << "Using r values: ";
    for (size_t i = 0, j = 0; i < ReadBatch::batches.size(); i++) {
      j = 0;
      for (auto r : ReadBatch::batches[i].rValues) {
        std::cerr << r << " (" << ReadBatch::batches[i].size << + ")";
        if ((i < ReadBatch::batches.size() - 1) || (j < ReadBatch::batches[i].rValues.size() - 1)) { std::cerr << ", "; }
        j++;
      }
    }
    std::cerr << '\n';
  }

  return true;
}

static Support testSequence(const Sequence& sequence) {
  int found = 0;
  int tests = 0;
  unsigned r = g_vanillaBloom->getKmerSize();
  if (sequence.size() >= r) {
    tests = sequence.size() - r + 1;

    int offset = 0;

    if (opt::error_correction) {
      stHashIterator it(sequence, g_spacedSeedsBloom->getSpacedSeeds(), SPACED_SEEDS_COUNT, SPACED_SEEDS_HASH_PER_SEED, r);
      for (; it != stHashIterator::end(); ++it, ++offset)
      {
        auto hitSeeds = g_spacedSeedsBloom->getHitSeeds(*it);
        if (hitSeeds.size() > 0) {
          it.snp({}, {}, g_vanillaBloom->getHashNum());
          if (g_vanillaBloom->contains(*it)) {
            found++;
          } else {
            bool success = false;
            for (auto hitSeed : hitSeeds) {
              for (auto seedIt = (hitSeed.begin() + std::round(hitSeed.size() * (1.00 - SPACED_SEEDS_SNP_FRACTION))); seedIt != hitSeed.end(); ++seedIt) {
                const auto pos = *seedIt;
                for (auto base : BASES) {
                  if (base == (unsigned char)(sequence[offset + pos])) { continue; }
                  it.snp({ pos }, { base }, g_vanillaBloom->getHashNum());
                  if (g_vanillaBloom->contains(*it)) {
                    success = true;
                    found++;
                    break;
                  }
                }
                if (success) { break; }
              }
              if (success) { break; }
            }
          }
        }
      }
    } else {
      for (ntHashIterator it(sequence, HASH_NUM, r); it != ntHashIterator::end(); ++it, ++offset)
      {
        if (g_vanillaBloom->contains(*it)) {
          found++;
        }
      }
    }
  }
  return Support(found, tests);
}

static Support testCombination(const std::string &head, const std::string &repeat, const std::string &tail, const int requestedTests)
{
  const auto windowSize = g_vanillaBloom->getKmerSize();

  auto plannedTests = requestedTests;
  if (plannedTests < opt::min_tests) { plannedTests = opt::min_tests; }
  if (plannedTests > opt::min_tests + MAX_TESTS_OFFSET) { return Support(); }

  const auto margin = getMargin(windowSize, plannedTests, repeat.size(), MIN_MARGIN);

  int possibleTests = head.size() + repeat.size() + tail.size() - windowSize + 1;
  if (possibleTests < plannedTests || long(head.size()) < margin || long(tail.size()) < margin) {
    return Support();
  }

  Sequence sequence;
  if (possibleTests > plannedTests + 1) {
    assert(long(head.size()) > margin || long(tail.size()) > margin);
    sequence = head.substr(head.size() - margin) + repeat + tail.substr(0, margin);
  } else {
    sequence = head + repeat + tail;
  }
  possibleTests = sequence.size() - windowSize + 1;

  assert(plannedTests <= possibleTests);
  assert(possibleTests <= plannedTests + 1);
  assert(int(sequence.size()) >= MIN_MARGIN + int(repeat.size()) + MIN_MARGIN);
  assert(int(sequence.size()) < int(windowSize) * 2);

  return testSequence(sequence);
}

static double expectedSpacingBetweenReads(const ContigPath& path) {
  const long pathLength = 100000;
  const double leftContigBaseCoverage = getContigBaseCoverage(path[0]);
  const double rightContigBaseCoverage = getContigBaseCoverage(path[path.size() - 1]);
  const double pathBaseCoverage = std::min(leftContigBaseCoverage, rightContigBaseCoverage);
  const double pathBases = pathBaseCoverage * (pathLength - opt::k + 1);

  double meanReadKmerContribution = 0;
  for (const auto& batch : ReadBatch::batches) {
    meanReadKmerContribution += batch.getFractionOfTotal() * (batch.size - opt::k + 1);
  }
  const double baseContributionRatio = ReadBatch::current.getFractionOfTotal() * (ReadBatch::current.size - opt::k + 1) / meanReadKmerContribution;

  const double approxNumOfReads = double(pathBases * baseContributionRatio) / double(opt::k * (ReadBatch::current.size - opt::k + 1));
  assert(approxNumOfReads > 2);

  return double(pathLength - ReadBatch::current.size) / double(approxNumOfReads - 1);
}

static Support determinePathSupport(const ContigPath& path)
{
  assert(path.size() >= 3);
  const Sequence repeat = getPathSequence(ContigPath(path.begin() + 1, path.end() - 1));
  const int repeatSize = repeat.size();
  assert(repeatSize >= 2);

  const long calculatedTests = std::round(expectedSpacingBetweenReads(path) * COV_APPROX_FORMULA_FACTOR);
  assert(calculatedTests >= 0);

  long requiredTests = calculatedTests;
  if (requiredTests < opt::min_tests) { requiredTests = opt::min_tests; }

  const int windowSize = g_vanillaBloom->getKmerSize();
  assert(windowSize >= 4);

  if (!windowLongEnough(windowSize, requiredTests, repeatSize, MIN_MARGIN)) {
    return Support(calculatedTests);
  }

  const auto &leftContig = path[0];
  const auto &rightContig = path[path.size() - 1];
  assert(windowSize >= MIN_MARGIN + repeatSize + MIN_MARGIN);

  const int leftDistance = distanceBetween(leftContig, path[1]);
  const int rightDistance = distanceBetween(path[path.size() - 2], rightContig);

  const int margin = getMargin(windowSize, requiredTests, repeat.size(), MIN_MARGIN);

  const auto heads = getTreeSequences(leftContig, -leftDistance, margin, false, opt::branching);
  const auto tails = getTreeSequences(rightContig, -rightDistance, margin, true, opt::branching);
  assert(heads.size() > 0);
  assert(tails.size() > 0);

  Support maxSupport(calculatedTests);
  bool unknown = false;

  const long combinations = heads.size() * tails.size();
  if (combinations > opt::branching * opt::branching) {
    unknown = true;
  } else {
    if (combinations >= PATH_COMBINATIONS_MULTITHREAD_THRESHOLD) {
      bool end = false;
      for (const auto head : heads) {
        #pragma omp critical (maxSupport)
        {
          if (unknown) { end = true; }
        }
        if (end) { break; }

        #pragma omp task firstprivate(head) shared(maxSupport, unknown)
        {
          bool end = false;
          for (const auto& tail : tails) {
            #pragma omp critical (maxSupport)
            {
              if (unknown) { end = true; }
            }
            if (end) { break; }

            auto support = testCombination(head, repeat, tail, requiredTests);

            #pragma omp critical (maxSupport)
            {
              if (support.unknown()) {
                unknown = true;
                end = true;
              } else if (support > maxSupport) {
                maxSupport = support;
              } else if (maxSupport.found == 0 && support.tests > maxSupport.tests) {
                maxSupport.tests = support.tests;
              }
            }
            if (end) { break; }
          }
        }
      }
    } else {
      for (const auto& head : heads) {
        if (unknown) { break; }

        for (const auto& tail : tails) {
          if (unknown) { break; }

          auto support = testCombination(head, repeat, tail, requiredTests);

          if (support.unknown()) {
            unknown = true;
            break;
          } else if (support > maxSupport) {
            maxSupport = support;
          } else if (maxSupport.found == 0 && support.tests > maxSupport.tests) {
            maxSupport.tests = support.tests;
          }
        }
      }
    }
  }

  if (combinations >= PATH_COMBINATIONS_MULTITHREAD_THRESHOLD) {
    #pragma omp taskwait
  }

  if (unknown) {
    return Support(calculatedTests);
  }

  maxSupport.calculatedTests = calculatedTests;
  return maxSupport;
}

static SupportMap
buildRepeatSupportMap(const ContigNode &repeat)
{
  SupportMap supportMap;
  bool unknown = false;
  in_edge_iterator inIt, inLast;
  for (std::tie(inIt, inLast) = in_edges(repeat, g_contigGraph);
       inIt != inLast; ++inIt)
  {
    const auto intig = source(*inIt, g_contigGraph);
    out_edge_iterator outIt, outLast;
    for (std::tie(outIt, outLast) = out_edges(repeat, g_contigGraph);
         outIt != outLast; ++outIt)
    {
      const auto outig = target(*outIt, g_contigGraph);
      auto support = determinePathSupport({ intig, repeat, outig });

      supportMap[intig.index()][outig.index()] = support;
      if (support.unknown()) {
        unknown = true;
      }
    }
  }

  if (unknown) {
    for (std::tie(inIt, inLast) = in_edges(repeat, g_contigGraph);
        inIt != inLast; ++inIt)
    {
      const auto intig = source(*inIt, g_contigGraph);
      out_edge_iterator outIt, outLast;
      for (std::tie(outIt, outLast) = out_edges(repeat, g_contigGraph);
          outIt != outLast; ++outIt)
      {
        const auto outig = target(*outIt, g_contigGraph);
        supportMap[intig.index()][outig.index()].reset();
      }
    }
  }

  return supportMap;
}

static void updateStats(Resolution& resolution, long &pathsKnown, long& pathsUnknown,
                        const SupportMap& repeatSupportMap,
                        bool inHistSample)
{
  for (const auto &intigIdxAndOutigsSupp : repeatSupportMap)
  {
    for (const auto &outigIdxAndSupp : intigIdxAndOutigsSupp.second)
    {
      const auto &support = outigIdxAndSupp.second;
      
      if (support.unknown()) {
        pathsUnknown++;
      } else {
        assert(support.found >= 0);
        assert(support.tests >= 0);
        pathsKnown++;
        if (inHistSample) {
          resolution.findsHistogram.insert(support.found);
          resolution.fractionFindsHistogram.insert(double(support.found) / double(support.tests));
        }
      }

      assert(support.calculatedTests >= 0);
      if (inHistSample) {
        resolution.calculatedTestsHistogram.insert(support.calculatedTests);
      }
    }
  }
}

static bool isSmallRepeat(const ContigNode& node)
{
  unsigned r = g_vanillaBloom->getKmerSize();
  return (!get(vertex_removed, g_contigGraph, node) && !node.sense() &&
      windowLongEnough(r, opt::min_tests, getContigSize(node), MIN_MARGIN) &&
      (in_degree(node, g_contigGraph) > 0 && out_degree(node, g_contigGraph) > 0) &&
      (in_degree(node, g_contigGraph) > 1 || out_degree(node, g_contigGraph) > 1));
}

static Resolution resolveRepeats()
{
  long total = (num_vertices(g_contigGraph) - num_vertices_removed(g_contigGraph)) / 2;
  long repeats = 0, pathsKnown = 0, pathsUnknown = 0;
  long pathsSupported = 0, pathsUnsupported = 0;

  progressStart("Path resolution (r = " + 
    std::to_string(g_vanillaBloom->getKmerSize()) + 
    ")", total * 2);

  Resolution resolution(ReadBatch::current, g_vanillaBloom->getKmerSize());

  Graph::vertex_iterator vertexStart, vertexEnd;
  boost::tie(vertexStart, vertexEnd) = vertices(g_contigGraph);

  iteratorMultithreading(vertexStart, vertexEnd,
  [&] (const ContigNode &node) {
    if (!get(vertex_removed, g_contigGraph, node)) {
      if (isSmallRepeat(node)) {
        return true;
      } else {
        #pragma omp critical (cerr)
        progressUpdate();
        return false;
      }
    }
    return false;
  },
  [&] (const ContigNode &node) {
    bool inHistSample;
    bool skip = false;
    #pragma omp critical(resolution)
    {
      repeats++;
      inHistSample = (repeats <= HIST_SAMPLE_SIZE);
      skip = (repeats > REPEAT_CASES_LIMIT);
    }

    if (!skip) {
      auto supportMap = buildRepeatSupportMap(node);

      #pragma omp critical(resolution)
      {
        resolution.repeatSupportMap[node.index()] = supportMap;
        updateStats(resolution, pathsKnown, pathsUnknown,
                    supportMap, inHistSample);
      }
    }

    #pragma omp critical (cerr)
    progressUpdate();
  });

  if (repeats > 0 && pathsKnown > 0) {
    for (const auto &findsAndCount : resolution.findsHistogram) {
      const auto &finds = findsAndCount.first;
      const auto &count = findsAndCount.second;
      if (finds >= opt::threshold) {
        pathsSupported += count;
      } else {
        pathsUnsupported += count;
      }
    }

    double sampleFactor = double(pathsKnown) / double(pathsSupported + pathsUnsupported);
    pathsSupported *= sampleFactor;
    pathsUnsupported *= sampleFactor;

    if (opt::verbose) {
      std::cerr << std::fixed;
      std::cerr << "Small repeats = " << repeats << "/" << total
        << " (" << double(repeats) / total * 100.0 << "%)\n";
      std::cerr << "Known support paths = " << pathsKnown << " / " << (pathsKnown + pathsUnknown)
        << " (" << double(pathsKnown) / (pathsKnown + pathsUnknown) * 100.0 << "%)\n";
      std::cerr << "Unknown support paths = " << pathsUnknown << " / " << (pathsKnown + pathsUnknown)
        << " (" << double(pathsUnknown) / (pathsKnown + pathsUnknown) * 100.0 << "%)\n";
      std::cerr << "Supported paths ~= " << pathsSupported << "/" << pathsKnown
        << " (" << double(pathsSupported) / pathsKnown * 100.0 << "%)\n";
      std::cerr << "Unsupported paths ~= " << pathsUnsupported << "/" << pathsKnown
        << " (" << double(pathsUnsupported) / pathsKnown * 100.0 << "%)\n";
      std::cerr << std::defaultfloat << std::flush;
    }

    if (double(pathsSupported) / double(pathsKnown) < SUPPORTED_PATHS_MIN) {
      std::cerr << "Insufficient support found. Is something wrong with the data?\n";
      resolution.failed = true;
    }
  } else {
    std::cerr << "No small junctions were found!\n";
    resolution.failed = true;
  }

  return resolution;
}

struct OldEdge {

  OldEdge(ContigNode u, ContigNode v):
    u(u), v(v) {} 

  ContigNode u, v;
};

struct NewEdge {

  NewEdge(ContigNode u, ContigNode v, Distance distance):
    u(u), v(v), distance(distance) {}

  ContigNode u, v;
  Distance distance;
};

struct NewVertex {

  NewVertex(ContigNode original, ContigNode node):
    original(original), node(node) {}

  ContigNode original, node;
};

class RepeatInstance {

public:

  RepeatInstance(const ContigNode instance,
                 const ContigNode original,
                 const std::vector<ContigNode> originalIntigs,
                 const std::vector<ContigNode> originalOutigs):
    instance(instance),
    original(original),
    originalIntigs(originalIntigs),
    originalOutigs(originalOutigs) {}

  bool inOriginalIntigs(const ContigNode &node) const {
    return std::find(originalIntigs.begin(), originalIntigs.end(), node)
            != originalIntigs.end();
  }

  bool inOriginalOutigs(const ContigNode &node) const {
    return std::find(originalOutigs.begin(), originalOutigs.end(), node)
            != originalOutigs.end();
  }

  RepeatInstance getReverse() const {
    std::vector<ContigNode> originalIntigsReverse;
    for (auto originalOutig : originalOutigs) {
      originalIntigsReverse.push_back(originalOutig ^ true);
    }

    std::vector<ContigNode> originalOutigsReverse;
    for (auto originalIntig : originalIntigs) {
      originalOutigsReverse.push_back(originalIntig ^ true);
    }

    return RepeatInstance(instance ^ true, original ^ true,
                          originalIntigsReverse, originalOutigsReverse);
  }

  const ContigNode instance;
  const ContigNode original;
  std::vector<ContigNode> originalIntigs;
  std::vector<ContigNode> originalOutigs;
  std::vector<std::reference_wrapper<const RepeatInstance>> intigsInstances;
  std::vector<std::reference_wrapper<const RepeatInstance>> outigsInstances;

};

static void processGraph(const Resolution &resolution,
                         ImaginaryContigPaths& supportedPaths,
                         ImaginaryContigPaths& unsupportedPaths)
{
  progressStart("New paths and vertices setup", resolution.repeatSupportMap.size() * 3);

  assert(!resolution.failed);

  std::vector<OldEdge> edges2remove;
  std::vector<NewEdge> edges2add;
  std::vector<NewVertex> vertices2add;

  std::map<int, std::vector<RepeatInstance>> repeatInstancesMap;

  size_t lastId = num_vertices(g_contigGraph) / 2;

  const auto start = resolution.repeatSupportMap.begin();
  const auto end = resolution.repeatSupportMap.end();

#ifdef _OPENMP
  int threads = omp_get_num_threads();
#else
  int threads = 1;
#endif

  // 1
  iteratorMultithreading(start, end,
  [&] (const std::pair<int, SupportMap> &repeatSupport) {
    (void)repeatSupport;
    return true;
  },
  [&] (const std::pair<int, SupportMap> &repeatSupport) {
    const auto repeat = ContigNode(repeatSupport.first);
    const auto &supportMap = repeatSupport.second;

    #pragma omp critical (repeatInstancesMap)
    {
      repeatInstancesMap.emplace(repeat.index(), std::vector<RepeatInstance>());
      repeatInstancesMap.emplace((repeat ^ true).index(), std::vector<RepeatInstance>());
    }

    for (const auto &intigIdxAndOutigsSupp : supportMap) {
      const auto intig = ContigNode(intigIdxAndOutigsSupp.first);
      for (const auto &outigIdxAndSupp : intigIdxAndOutigsSupp.second) {
        const auto outig = ContigNode(outigIdxAndSupp.first);
        const auto &support = outigIdxAndSupp.second;
        int dist1 = distanceBetween(intig, repeat);
        int dist2 = distanceBetween(repeat, outig);

        ImaginaryContigPath path = { { intig, 0 }, { repeat, dist1 }, { outig, dist2 } };

        if (support.good()) {
          #pragma omp critical (supportedPaths)
          supportedPaths.insert(path);
        } else {
          #pragma omp critical (unsupportedPaths)
          unsupportedPaths.insert(path);

          #pragma omp critical (supportedPaths)
          {
            if (supportedPaths.find(path) != supportedPaths.end()) {
              supportedPaths.erase(path);
            }
          }
        }
      }
    }

    #pragma omp critical (cerr)
    progressUpdate();
  }, std::min(4, threads));

  // 2
  for (auto it = start; it != end; it++) {
    const std::pair<int, SupportMap> &repeatSupport = *it;
    const auto repeat = ContigNode(repeatSupport.first);
    const auto &supportMap = repeatSupport.second;

    assert(repeatInstancesMap.find(repeat.index()) != repeatInstancesMap.end());
    assert(repeatInstancesMap.find((repeat ^ true).index()) != repeatInstancesMap.end());

    auto &repeatInstances = repeatInstancesMap.at(repeat.index());
    auto &repeatInstancesReverse = repeatInstancesMap.at((repeat ^ true).index());

    assert(repeatInstances.size() == 0);
    assert(repeatInstancesReverse.size() == 0);

    for (const auto &intigIdxAndOutigsSupp : supportMap) {
      const auto intig = ContigNode(intigIdxAndOutigsSupp.first);
      const auto &outigsSupp = intigIdxAndOutigsSupp.second;

      std::vector<ContigNode> supportedOutigs;
      for(const auto &outigIdxAndSupp : outigsSupp) {
        const auto &outig = ContigNode(outigIdxAndSupp.first);
        const auto &support = outigIdxAndSupp.second;
        if (support.good()) {
          supportedOutigs.push_back(outig);
        }
      }

      bool matched = false;
      for (auto &instance : repeatInstances) {
        if (instance.originalOutigs.size() == supportedOutigs.size()) {
          matched = true;
          for (const auto &outig: supportedOutigs) {
            bool found = false;
            for (const auto &instanceOutig : instance.originalOutigs) {
              if (outig == instanceOutig) {
                found = true;
                break;
              }
            }
            if (!found) {
              matched = false;
              break;
            }
          }
        }
        if (matched) {
          instance.originalIntigs.push_back(intig);
          break;
        }
      }

      if (!matched) {
        if (supportedOutigs.size() > 0) {
          std::vector<ContigNode> intigs = { intig };
          if (repeatInstances.size() == 0) {
            repeatInstances.push_back(RepeatInstance(repeat, repeat, intigs, supportedOutigs));
          } else {
            ContigNode repeatCopy = ContigNode(lastId++, repeat.sense());
            repeatInstances.push_back(RepeatInstance(repeatCopy, repeat, intigs, supportedOutigs));
          }
        }
      }
    }

    if (repeatInstances.size() > 0) {
      std::set<int> intigIdxs;

      for (const auto &instance : repeatInstances) {
        for (const auto &intig : instance.originalIntigs) {
          assert(intigIdxs.find(intig.index()) == intigIdxs.end());
          intigIdxs.insert(intig.index());
        }
        assert(instance.originalOutigs.size() > 0);
        repeatInstancesReverse.push_back(instance.getReverse());
      }
    } else {
      auto instance = RepeatInstance(repeat, repeat, {}, {});
      repeatInstances.push_back(instance);
      repeatInstancesReverse.push_back(instance.getReverse());
    }

    progressUpdate();
  }

  // 3
  iteratorMultithreading(start, end,
  [&] (const std::pair<int, SupportMap> &repeatSupport) {
    (void)repeatSupport;
    return true;
  },
  [&] (const std::pair<int, SupportMap> &repeatSupport) {
    const auto repeat = ContigNode(repeatSupport.first);

    auto &repeatInstances = repeatInstancesMap.at(repeat.index());

    std::list<RepeatInstance> tempInstances;

    for (auto &instance : repeatInstances) {
      for (const auto &intig : instance.originalIntigs) {
        if (repeatInstancesMap.find(intig.index()) != repeatInstancesMap.end()) {
          const auto &intigInstances = repeatInstancesMap.at(intig.index());
          for (const auto &intigInstance : intigInstances) {
            if (intigInstance.inOriginalOutigs(repeat)) {
              instance.intigsInstances.push_back(intigInstance);
            }
          }
        } else {
          tempInstances.push_back(RepeatInstance(intig, intig, {}, {}));
          instance.intigsInstances.push_back(tempInstances.back());
        }
      }

      for (const auto &outig : instance.originalOutigs) {
        if (repeatInstancesMap.find(outig.index()) != repeatInstancesMap.end()) {
          const auto &outigInstances = repeatInstancesMap.at(outig.index());
          for (const auto &outigInstance : outigInstances) {
            if (outigInstance.inOriginalIntigs(repeat)) {
              instance.outigsInstances.push_back(outigInstance);
            }
          }
        } else {
          tempInstances.push_back(RepeatInstance(outig, outig, {}, {}));
          instance.outigsInstances.push_back(tempInstances.back());
        }
      }

      if (instance.instance == instance.original) {
        in_edge_iterator inIt, inLast;
        for (std::tie(inIt, inLast) = in_edges(instance.original, g_contigGraph);
          inIt != inLast; ++inIt)
        {
          #pragma omp critical (edges2remove)
          edges2remove.push_back(OldEdge(source(*inIt, g_contigGraph), instance.original));
        }

        out_edge_iterator outIt, outLast;
        for (std::tie(outIt, outLast) = out_edges(instance.original, g_contigGraph);
          outIt != outLast; ++outIt)
        {
          #pragma omp critical (edges2remove)
          edges2remove.push_back(OldEdge(instance.original, target(*outIt, g_contigGraph)));
        }
      } else {
        #pragma omp critical (vertices2add)
        vertices2add.push_back(NewVertex(instance.original, instance.instance));
      }
      
      for (const RepeatInstance &intigInstance : instance.intigsInstances) {
        #pragma omp critical (edges2add)
        edges2add.push_back(NewEdge(intigInstance.instance, instance.instance,
          get(edge_bundle, g_contigGraph, edge(intigInstance.original, instance.original, g_contigGraph).first)));
      }

      for (const RepeatInstance &outigInstance : instance.outigsInstances) {
        #pragma omp critical (edges2add)
        edges2add.push_back(NewEdge(instance.instance, outigInstance.instance,
          get(edge_bundle, g_contigGraph, edge(instance.original, outigInstance.original, g_contigGraph).first)));
      }
    }

    #pragma omp critical (cerr)
    progressUpdate();
  });

  std::sort(vertices2add.begin(), vertices2add.end(),
  [] (const NewVertex& v1, const NewVertex& v2) {
    return v1.node.index() < v2.node.index();
  });
  std::sort(edges2add.begin(), edges2add.end(),
  [] (const NewEdge &e1, const NewEdge &e2) {
    return (e1.u.index() < e2.u.index()) || (e1.u.index() == e2.u.index() && e1.v.index() < e2.v.index());
  });

  int modifications = edges2remove.size() + vertices2add.size() + edges2add.size();
  progressStart("Graph modification", modifications);

  g_contigNames.unlock();
  for (const auto &oldEdge : edges2remove) {
    if (g_contigGraph.edge(oldEdge.u, oldEdge.v).second) {
      remove_edge(oldEdge.u, oldEdge.v, g_contigGraph);
    }

    progressUpdate();
  }
  for (const auto &newVertex : vertices2add) {
    assert(in_degree(newVertex.original, g_contigGraph) == 0);
    assert(out_degree(newVertex.original, g_contigGraph) == 0);

    assert(g_contigSequences.size() == newVertex.node.index());
    assert(g_contigComments.size() == newVertex.node.id());

    g_contigSequences.push_back(getContigSequence(newVertex.original));
    g_contigSequences.push_back(getContigSequence(newVertex.original ^ true));

    std::string name = createContigName();
    put(vertex_name, g_contigGraph, newVertex.node, name);
    add_vertex(get(vertex_bundle, g_contigGraph, newVertex.original), g_contigGraph);

    g_contigComments.push_back(getContigComment(newVertex.original));

    assert(in_degree(newVertex.node, g_contigGraph) == 0);
    assert(out_degree(newVertex.node, g_contigGraph) == 0);

    progressUpdate();
  }
  for (const auto &newEdge : edges2add) {
    if (!g_contigGraph.edge(newEdge.u, newEdge.v).second) {
      g_contigGraph.add_edge(newEdge.u, newEdge.v, newEdge.distance);
    }

    progressUpdate();
  }
  g_contigNames.lock();
}

void writeHistograms(const Resolution& resolution, const std::string& prefix, int subiteration) {
  if (opt::verbose) {
    std::cerr << "Writing algorithm histograms..." << std::flush;
  }

  std::string findsFilename = prefix + "-r" + std::to_string(resolution.r) + "-" + std::to_string(subiteration + 1) + "-finds.tsv";
  std::ofstream findsFile(findsFilename.c_str());
  findsFile << resolution.findsHistogram;

  std::string fractionFindsFilename = prefix + "-r" + std::to_string(resolution.r) + "-" + std::to_string(subiteration + 1) + "-percent-finds.tsv";
  std::ofstream fractionFindsFile(fractionFindsFilename.c_str());
  fractionFindsFile << resolution.fractionFindsHistogram;

  std::string calculatedTestsFilename = prefix + "-r" + std::to_string(resolution.r) + "-" + std::to_string(subiteration + 1) + "-calculated-tests.tsv";
  std::ofstream calculatedTestsFile(calculatedTestsFilename.c_str());
  calculatedTestsFile << resolution.calculatedTestsHistogram;

  if (opt::verbose) {
    std::cerr << " Done!" << std::endl;
  }
}

void
resolveShort(const std::vector<std::string>& readFilepaths,
             ImaginaryContigPaths& supportedPaths,
             ImaginaryContigPaths& unsupportedPaths)
{
  if (!determineShortReadStats(readFilepaths)) {
    return;
  }

  if (opt::verbose) { std::cerr << "\nRunning resolution algorithm...\n"; }

  assert(g_contigSequences.size() > 0);
  assert(g_contigSequences.size() / 2 == g_contigComments.size());
  assert(ReadBatch::batches.size() > 0);

  std::vector<std::pair<int, Histogram>> histograms;
  for (auto batch : ReadBatch::batches) {
    ReadBatch::current = batch;

    for (int r : ReadBatch::current.rValues) {
      if (int(r) < int(opt::k)) {
        std::cerr << "r value " << r << "(" << ReadBatch::current.size << ") is too short - skipping." << std::endl;
        continue;
      }

      if (opt::verbose) { std::cerr << "\nRead size = " << batch.size << ", r = " << r << " ...\n\n"; }

      buildFilters(readFilepaths, r, opt::bloomSize);

      for (size_t j = 0; j < MAX_SUBITERATIONS; j++) {
        if (opt::verbose) { std::cerr << "\nSubiteration " << j + 1 << "...\n"; }

        int unsupportedCountPrev = unsupportedPaths.size();

        Resolution resolution = resolveRepeats();

        if (!resolution.failed) {
          processGraph(resolution, supportedPaths, unsupportedPaths);
          assembleContigs();
          if (!opt::histPrefix.empty()) {
            writeHistograms(resolution, opt::histPrefix, j);
          }
        }

        int newUnsupportedCount = unsupportedPaths.size() - unsupportedCountPrev;
        assert(newUnsupportedCount >= 0);
        if (newUnsupportedCount == 0) { break; }
      }
    }
  }

  if (opt::verbose) { std::cerr << "Resolution algorithm done.\n\n"; }
}
