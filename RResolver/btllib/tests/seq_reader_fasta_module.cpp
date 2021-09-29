#include "btllib/seq_reader.hpp"
#include "helpers.hpp"

#include <iostream>

int
main()
{
  const char* ids[] = { "asdf", "ghjk" };
  const char* seqs[] = { "ACTG", "TGCA" };

  for (int iteration = 0; iteration < 3; iteration++) {
    std::cerr << "Iteration " << iteration + 1 << std::endl;

    std::cerr << "Test small FASTA" << std::endl;
    btllib::SeqReader reader("../tests/input.fa.gz.bz2.xz.lrz",
                             btllib::SeqReader::Flag::SHORT_MODE);
    TEST_ASSERT_EQ(reader.get_format(), btllib::SeqReader::Format::FASTA)

    size_t i = 0;
    for (const auto record : reader) {
      TEST_ASSERT_EQ(record.id, ids[i]);
      TEST_ASSERT_EQ(record.seq, seqs[i])
      i++;
    }
  }

  return 0;
}