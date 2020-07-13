#ifndef BTLLIB_TESTS_HELPERS
#define BTLLIB_TESTS_HELPERS

#include <random>
#include <chrono>
#include <string>
#include <functional>

inline std::string
get_random_sequence(size_t size) {
  static std::default_random_engine generator(std::chrono::system_clock::now().time_since_epoch().count() + 9999999);
  static std::uniform_int_distribution<int> distribution_actg(0, 3);
  static auto gen_random_actg = std::bind(distribution_actg, generator);
  std::string seq;
  for (size_t i = 0; i < size; i++) {
    seq += "ACTG"[(gen_random_actg())];
  }
  return seq;
}

inline std::string
get_random_name(size_t size) {
  static std::default_random_engine generator(std::chrono::system_clock::now().time_since_epoch().count() + 9999999);
  static std::uniform_int_distribution<int> distribution_alphabet(65, 90);
  static auto gen_random_alphabet = std::bind(distribution_alphabet, generator);
  std::string name;
  for (size_t i = 0; i < size; i++) {
    name += char(gen_random_alphabet());
  }
  return name;
}


#endif