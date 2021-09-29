#include "btllib/seq_reader.hpp"
#include "helpers.hpp"

#include <fstream>
#include <iostream>

int
main()
{
  // Test GFA2
  /*std::cerr << "Test GFA2" << std::endl;
  btllib::SeqReader reader_gfa2("../tests/input.gfa2",
                                btllib::SeqReader::Flag::SHORT_MODE |
                                  btllib::SeqReader::Flag::FOLD_CASE);
  TEST_ASSERT_EQ(reader_gfa2.get_format(), btllib::SeqReader::Format::GFA2);

  i = 0;
  while ((record = reader_gfa2.read())) {
    TEST_ASSERT_EQ(record.seq, seqs[i]);
    TEST_ASSERT(record.qual.empty());

    i++;
  }
  TEST_ASSERT_EQ(i, 2);*/

  return 0;
}