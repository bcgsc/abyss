#ifndef RUTILS_H
#define RUTILS_H 1

#include <vector>
#include <string>
#include <algorithm>

#if _OPENMP
const unsigned MAX_SIMULTANEOUS_TASKS = 60000;
# include <omp.h>
#endif

const double PROGRESS_PRINT_FRACTION = 0.01;

template <typename IteratorT, typename FilterT, typename ActionT>
void iteratorMultithreading(const IteratorT& start, const IteratorT &end,
                        const FilterT &filter, const ActionT &action,
                        const int threads = 0)
{
#if _OPENMP
  int threadsBefore = 1;
  if (threads > 0) {
    threadsBefore = omp_get_max_threads();
    omp_set_num_threads(threads);
  }
#endif
  IteratorT it = start;
  while (it != end) {
    #pragma omp parallel
    #pragma omp single
    {
#if _OPENMP
      for (unsigned i = 0; i < MAX_SIMULTANEOUS_TASKS; i++)
#endif
      {
        if (filter(*it)) {
          #pragma omp task firstprivate (it)
          {
            action(*it);
          }
        }
        ++it;
#if _OPENMP
        if (it == end) break;
#endif
      }
    }
  }
#if _OPENMP
  if (threads > 0) {
    omp_set_num_threads(threadsBefore);
  }
#endif
}

void progressStart(const std::string& name, unsigned total);
void progressUpdate();

unsigned nChoosek(unsigned n, unsigned k);
std::vector<std::vector<unsigned>> genPermutations(unsigned n);
std::vector<std::vector<unsigned>> genCombinations(unsigned n, unsigned k);

#endif
