#include "btllib/seq_reader.hpp"
#include "helpers.hpp"

#include <fstream>
#include <iostream>

int
main()
{
  const char* ids[] = { "q1", "q2" };
  const char* seqs[] = { "ACTG", "TGCA" };
  const char* quals[] = { "!@^&", "(#&$" };

  for (int iteration = 0; iteration < 3; iteration++) {
    std::cerr << "Iteration " << iteration + 1 << std::endl;

    std::cerr << "Test small SAM" << std::endl;
    btllib::SeqReader reader("../tests/input.bam",
                             btllib::SeqReader::Flag::SHORT_MODE);
    TEST_ASSERT_EQ(reader.get_format(), btllib::SeqReader::Format::SAM);

    size_t i = 0;
    for (const auto record : reader) {
      TEST_ASSERT_EQ(record.id, ids[i])
      TEST_ASSERT_EQ(record.seq, seqs[i])
      TEST_ASSERT_EQ(record.qual, quals[i])
      i++;
    }

    std::cerr << "Test larger SAM file" << std::endl;
    btllib::SeqReader large_fastq_reader("../tests/large.fq",
                                         btllib::SeqReader::Flag::SHORT_MODE);
    btllib::SeqReader large_bam_reader("../tests/large.bam",
                                       btllib::SeqReader::Flag::SHORT_MODE);
    while (true) {
      btllib::SeqReader::Record record1 = large_fastq_reader.read();
      btllib::SeqReader::Record record2 = large_bam_reader.read();
      if (!bool(record1) && !bool(record2)) {
        break;
      }
      TEST_ASSERT(record1)
      TEST_ASSERT(record2)
      TEST_ASSERT_EQ(record1.id, record2.id);
      TEST_ASSERT_EQ(record1.comment, record2.comment);
      TEST_ASSERT_EQ(record1.seq, record2.seq);
      TEST_ASSERT_EQ(record1.qual, record2.qual);
    }
  }

  return 0;
}