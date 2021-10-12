#include "btllib/util.hpp"
#include "helpers.hpp"

#include <iostream>

int
main()
{
  std::string teststring("  AC  TG ");
  btllib::trim(teststring);
  TEST_ASSERT_EQ(teststring, "AC  TG");

  btllib::CString testcstring("  TG AC    ");
  btllib::trim(testcstring);
  TEST_ASSERT_EQ(std::string(testcstring), "TG AC");

  auto actg_split = btllib::split("A|C|T|G", "|");
  TEST_ASSERT_EQ(actg_split.size(), 4);
  TEST_ASSERT_EQ(actg_split[0], "A");
  TEST_ASSERT_EQ(actg_split[1], "C");
  TEST_ASSERT_EQ(actg_split[2], "T");
  TEST_ASSERT_EQ(actg_split[3], "G");

  auto actg_join = btllib::join(actg_split, "/");
  TEST_ASSERT_EQ(actg_join, "A/C/T/G");

  TEST_ASSERT(btllib::startswith(actg_join, "A/C"));
  TEST_ASSERT(btllib::endswith(actg_join, "/T/G"));

  return 0;
}