#include "btllib/seq.hpp"
#include "helpers.hpp"

#include <string>

int
main()
{
  std::string seq = "ACGTACACTGGACTGAGTCT";
  std::string rc = "AGACTCAGTCCAGTGTACGT";
  TEST_ASSERT_EQ(btllib::get_reverse_complement(seq), rc);
  return 0;
}