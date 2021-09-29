#ifndef BTLLIB_TESTS_HELPERS
#define BTLLIB_TESTS_HELPERS

#include <chrono>
#include <cstdlib>
#include <functional>
#include <random>
#include <string>

#define TEST_ASSERT(x)                                                         \
  if (!(x)) {                                                                  \
    std::cerr << __FILE__ ":" << __LINE__ << ":TEST_ASSERT: " #x " is false"   \
              << std::endl;                                                    \
    std::exit(EXIT_FAILURE);                                                   \
  }

#define TEST_ASSERT_RELATIONAL(x, y, op)                                       \
  if (!((x)op(y))) {                                                           \
    std::cerr << __FILE__ ":" << __LINE__                                      \
              << ":TEST_ASSERT_RELATIONAL: " #x " " #op " " #y << '\n'         \
              << #x " = " << x << '\n'                                         \
              << #y " = " << y << std::endl;                                   \
    std::exit(EXIT_FAILURE);                                                   \
  }

#define TEST_ASSERT_EQ(x, y) TEST_ASSERT_RELATIONAL(x, y, ==)
#define TEST_ASSERT_NE(x, y) TEST_ASSERT_RELATIONAL(x, y, !=)
#define TEST_ASSERT_GE(x, y) TEST_ASSERT_RELATIONAL(x, y, >=)
#define TEST_ASSERT_GT(x, y) TEST_ASSERT_RELATIONAL(x, y, >)
#define TEST_ASSERT_LE(x, y) TEST_ASSERT_RELATIONAL(x, y, <=)
#define TEST_ASSERT_LT(x, y) TEST_ASSERT_RELATIONAL(x, y, <)

inline int
get_random(int min, int max)
{
  static std::default_random_engine random_generator(
    std::chrono::system_clock::now().time_since_epoch().count() + 9999999);
  std::uniform_int_distribution<int> distribution(min, max);
  return distribution(random_generator);
}

inline std::string
get_random_seq(size_t size)
{
  static std::default_random_engine random_generator(
    std::chrono::system_clock::now().time_since_epoch().count() + 9999999);
  static std::uniform_int_distribution<int> distribution_actg(0, 3);
  static auto gen_random_actg = std::bind(distribution_actg, random_generator);
  std::string seq;
  for (size_t i = 0; i < size; i++) {
    seq += "ACTG"[gen_random_actg()];
  }
  return seq;
}

inline std::string
split_seq_multiline(std::string seq)
{
  static std::default_random_engine random_generator2(
    std::chrono::system_clock::now().time_since_epoch().count() + 9999998);
  static std::uniform_real_distribution<double> distribution_newline(0.0, 1.0);
  static auto gen_random_newline =
    std::bind(distribution_newline, random_generator2);
  for (size_t i = 1; i < seq.size(); i++) {
    if (gen_random_newline() <= 0.012 && seq[i - 1] != '\n') {
      seq.insert(i, "\n");
    }
  }
  return seq;
}

inline std::string
get_random_name(size_t size)
{
  static std::default_random_engine random_generator(
    std::chrono::system_clock::now().time_since_epoch().count() + 9999999);
  static std::uniform_int_distribution<int> distribution_alphabet(65, 90);
  static auto gen_random_alphabet =
    std::bind(distribution_alphabet, random_generator);
  std::string name;
  for (size_t i = 0; i < size; i++) {
    name += char(gen_random_alphabet());
  }
  return name;
}

#endif