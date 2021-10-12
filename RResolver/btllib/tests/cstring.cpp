#include "btllib/cstring.hpp"
#include "helpers.hpp"

#include <iostream>

int
main()
{
  btllib::CString cstring("ACTG");

  TEST_ASSERT_EQ(cstring[2], 'T');
  cstring[2] = 'C';
  TEST_ASSERT_EQ(cstring[2], 'C');

  cstring.erase(1, 2);
  TEST_ASSERT_EQ(std::string(cstring), "AG");

  cstring.erase();
  TEST_ASSERT_EQ(cstring.size(), 0);
  TEST_ASSERT_EQ(std::string(cstring), "");

  return 0;
}