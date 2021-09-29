#include "btllib/data_stream.hpp"

#include "helpers.hpp"

#include <cstdio>
#include <cstring>
#include <fstream>
#include <thread>

int
main()
{
  const char* txt = "data_stream test";
  char* line = new char[128];
  size_t line_len;

  // Test .gz
  const auto gz_filename = get_random_name(64) + ".gz";

  std::cerr << "Test .gz write" << std::endl;
  btllib::DataSink gz_sink(gz_filename, false);
  TEST_ASSERT_EQ(fwrite(txt, 1, strlen(txt), gz_sink), strlen(txt));
  gz_sink.close();

  std::cerr << "Test .gz read" << std::endl;
  btllib::DataSource gz_source(gz_filename);
  TEST_ASSERT_GT(getline(&line, &line_len, gz_source), 0);
  gz_source.close();
  TEST_ASSERT_EQ(strcmp(line, txt), 0);

  std::remove(gz_filename.c_str());

  // Test .xz
  const auto xz_filename = get_random_name(64) + ".xz";

  std::cerr << "Test .xz write" << std::endl;
  btllib::DataSink xz_sink(xz_filename, false);
  TEST_ASSERT_EQ(fwrite(txt, 1, strlen(txt), xz_sink), strlen(txt));
  xz_sink.close();

  std::cerr << "Test .xz read" << std::endl;
  btllib::DataSource xz_source(xz_filename);
  TEST_ASSERT_GT(getline(&line, &line_len, xz_source), 0);
  xz_source.close();
  TEST_ASSERT_EQ(strcmp(line, txt), 0);

  std::remove(xz_filename.c_str());

  // Test .lrz
  const auto lrz_filename = get_random_name(64) + ".lrz";

  std::cerr << "Test .lrz write" << std::endl;
  btllib::DataSink lrz_sink(lrz_filename, false);
  TEST_ASSERT_EQ(fwrite(txt, 1, strlen(txt), lrz_sink), strlen(txt));
  lrz_sink.close();

  std::cerr << "Test .lrz read" << std::endl;
  btllib::DataSource lrz_source(lrz_filename);
  TEST_ASSERT_GT(getline(&line, &line_len, lrz_source), 0);
  lrz_source.close();
  TEST_ASSERT_EQ(strcmp(line, txt), 0);

  std::remove(lrz_filename.c_str());

  return 0;
}