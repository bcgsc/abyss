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
  auto gz_filename = get_random_name(64) + ".gz";

  std::cerr << "Test .gz write" << std::endl;
  auto gz_sink = btllib::DataSink(gz_filename, false);
  fwrite(txt, strlen(txt), 1, gz_sink);
  gz_sink.close();

  std::cerr << "Test .gz read" << std::endl;
  auto gz_source = btllib::DataSource(gz_filename);
  getline(&line, &line_len, gz_source);
  gz_source.close();
  assert(strcmp(line, txt) == 0);

  std::remove(gz_filename.c_str());

  // Test .xz
  auto xz_filename = get_random_name(64) + ".xz";

  std::cerr << "Test .xz write" << std::endl;
  auto xz_sink = btllib::DataSink(xz_filename, false);
  fwrite(txt, strlen(txt), 1, xz_sink);
  xz_sink.close();

  std::cerr << "Test .xz read" << std::endl;
  auto xz_source = btllib::DataSource(xz_filename);
  getline(&line, &line_len, xz_source);
  xz_source.close();
  assert(strcmp(line, txt) == 0);

  std::remove(xz_filename.c_str());

  return 0;
}