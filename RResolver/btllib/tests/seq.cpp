#include "btllib/seq.hpp"

#include <cassert>
#include <string>

int
main()
{
  std::string seq = "ACGTACACTGGACTGAGTCT";
  std::string rc = "AGACTCAGTCCAGTGTACGT";
  assert(btllib::get_reverse_complement(seq) == rc);
  return 0;
}