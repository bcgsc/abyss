#include "RUtils.h"

#include "Common/Options.h"

#include <iostream>
#include <cassert>
#include <iomanip>

static std::string progressName;
static unsigned progressNumber = 0;
static unsigned progressTotal;
static unsigned progressUpdates;
static unsigned progressLastPrinted;
static bool progressRunning = false;

void progressStart(const std::string& name, unsigned total) {
  if (opt::verbose) {
    progressName = name;
    std::string tempName = name;
    tempName[0] = std::tolower(tempName[0]);
    std::cerr << '\n' << ++progressNumber << ". Starting " + tempName << "...\n";
    std::cerr << "Progress: 0%" << std::flush;
    progressTotal = total;
    progressUpdates = 0;
    progressLastPrinted = 0;
    assert(!progressRunning);
    progressRunning = true;
  }
}

void progressUpdate() {
  if (opt::verbose) {
    assert(progressRunning);
    progressUpdates++;
    assert(progressUpdates <= progressTotal);
    double fraction = double(progressUpdates) / progressTotal;

    if (progressUpdates == progressTotal) {
      std::cerr << "\rProgress: 100%\n" << progressName << " done." << std::endl;
      progressRunning = false;
    } else if (double(progressUpdates - progressLastPrinted) / progressTotal 
                >= PROGRESS_PRINT_FRACTION && progressUpdates < progressTotal)
    {
      std::cerr << "\rProgress: " << int(fraction * 100.0) << "%" << std::flush;
      progressLastPrinted = progressUpdates;
    }
  }
}

unsigned nChoosek(unsigned n, unsigned k) {
  if (k > n) return 0;
  if (k * 2 > n) k = n - k;
  if (k == 0) return 1;

  long result = n;
  for(unsigned i = 2; i <= k; ++i) {
    result *= (n - i + 1);
    result /= i;
  }
  return result;
}

std::vector<std::vector<unsigned>> genPermutations(unsigned n) {
  std::vector<std::vector<unsigned>> permutations;
  std::vector<unsigned> permutation;
  for (unsigned i = 0; i < n; i++) {
    permutation.push_back(i);
  }
  for (unsigned i = 0; i < n; i++) {
    permutations.push_back(permutation);
    std::next_permutation(permutation.begin(), permutation.end());
  }
  return permutations;
}

std::vector<std::vector<unsigned>> genCombinations(unsigned n, unsigned k) {
  std::vector<bool> v(n);
  std::fill(v.begin(), v.begin() + k, true);
  std::vector<std::vector<unsigned>> combinations;
  do {
    std::vector<unsigned> combination;
    for (unsigned i = 0; i < n; ++i) {
      if (v[i]) {
        combination.push_back(i);
      }
    }
    combinations.push_back(combination);
  } while (std::prev_permutation(v.begin(), v.end()));
  return combinations;
}
