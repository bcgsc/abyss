#include "btllib/seq_reader.hpp"
#include "helpers.hpp"

#include <fstream>
#include <iostream>

int
main()
{
  const char* ids[] = { "asdf", "ghjk" };
  const char* seqs[] = { "ACTG", "TGCA" };
  std::string random_filename;

  for (int iteration = 0; iteration < 3; iteration++) {
    std::cerr << "Iteration " << iteration + 1 << std::endl;

    std::cerr << "Test multiline FASTA" << std::endl;
    btllib::SeqReader reader("../tests/input_multiline.fa",
                             btllib::SeqReader::Flag::SHORT_MODE);

    TEST_ASSERT_EQ(reader.get_format(), btllib::SeqReader::Format::FASTA);

    size_t i = 0;
    for (const auto record : reader) {
      TEST_ASSERT_EQ(record.id, ids[i]);
      TEST_ASSERT_EQ(record.seq, seqs[i]);
      i++;
    }
    TEST_ASSERT_EQ(i, 2);
    reader.close();

    std::cerr << "Test random multiline FASTA file" << std::endl;
    std::vector<std::string> generated_ids;
    std::vector<std::string> generated_comments;
    std::vector<std::string> generated_seqs;
    random_filename = get_random_name(64);
    std::ofstream random_seqs_multiline(random_filename);
    for (int s = 0; s < 500; s++) {
      std::string id, comment_spaces, comment, seq, qual;

      id = get_random_name(10);
      comment_spaces = std::string(get_random(1, 10), ' ');
      comment = get_random_name(20);
      size_t seq_size = get_random(100, 2000);
      seq = get_random_seq(seq_size);

      std::string seq_multiline = split_seq_multiline(seq);
      std::string newline = get_random(0, 1) == 1 ? "\r\n" : "\n";

      random_seqs_multiline << '>' << id << ' ' << comment << newline
                            << seq_multiline << newline;

      generated_ids.push_back(id);
      generated_comments.push_back(comment);
      generated_seqs.push_back(seq);
    }
    random_seqs_multiline.close();
    std::cerr << random_filename << std::endl;

    btllib::SeqReader random_reader(random_filename,
                                    btllib::SeqReader::Flag::SHORT_MODE);
    TEST_ASSERT_EQ(random_reader.get_format(),
                   btllib::SeqReader::Format::FASTA);
    btllib::SeqReader::Record record;
    for (i = 0; (record = random_reader.read()); i++) {
      TEST_ASSERT_EQ(record.id, generated_ids[i]);
      TEST_ASSERT_EQ(record.comment, generated_comments[i]);
      TEST_ASSERT_EQ(record.seq, generated_seqs[i]);
    }
    TEST_ASSERT_EQ(i, 500);

    random_reader.close();
    std::remove(random_filename.c_str());
  }

  return 0;
}