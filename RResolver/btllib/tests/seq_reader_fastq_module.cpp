#include "btllib/seq_reader.hpp"
#include "helpers.hpp"

#include <fstream>
#include <iostream>

int
main()
{
  const char* ids[] = { "asdf", "ghjk" };
  const char* seqs[] = { "ACTG", "TGCA" };
  const char* quals[] = { "!@^&", "(#&$" };
  std::string random_filename;

  for (int iteration = 0; iteration < 3; iteration++) {
    std::cerr << "Iteration " << iteration + 1 << std::endl;

    size_t i;
    btllib::SeqReader::Record record;

    std::cerr << "Test small FASTA and FASTQ simultaneously" << std::endl;
    btllib::SeqReader reader_fasta("../tests/input.fa.gz.bz2.xz.lrz",
                                   btllib::SeqReader::Flag::SHORT_MODE);
    btllib::SeqReader reader_fastq("../tests/input.fq.tar.xz",
                                   btllib::SeqReader::Flag::SHORT_MODE);

    TEST_ASSERT_EQ(reader_fasta.get_format(), btllib::SeqReader::Format::FASTA);
    TEST_ASSERT_EQ(reader_fastq.get_format(), btllib::SeqReader::Format::FASTQ);

    i = 0;
    bool success_fasta = false, success_fastq = false;
    for (;;) {
      if ((success_fasta = (record = reader_fasta.read()))) {
        TEST_ASSERT_EQ(record.id, ids[i]);
        TEST_ASSERT_EQ(record.seq, seqs[i]);
        TEST_ASSERT(record.qual.empty());
      }

      if ((success_fastq = (record = reader_fastq.read()))) {
        TEST_ASSERT_EQ(record.id, ids[i]);
        TEST_ASSERT_EQ(record.seq, seqs[i]);
        TEST_ASSERT_EQ(record.qual, quals[i]);
      }

      if (++i == 3) {
        TEST_ASSERT(success_fasta == false && success_fastq == false);
        break;
      } else {
        TEST_ASSERT(success_fasta == true && success_fastq == true);
      }
    }

    std::cerr << "Test larger randomly generated FASTQ file" << std::endl;
    std::vector<std::string> generated_ids;
    std::vector<std::string> generated_comments;
    std::vector<std::string> generated_seqs;
    std::vector<std::string> generated_quals;
    random_filename = get_random_name(64);
    std::ofstream random_seqs(random_filename);
    for (int s = 0; s < 500; s++) {
      std::string id, comment_spaces, comment, seq, qual;

      id = get_random_name(10);
      comment_spaces = std::string(get_random(1, 10), ' ');
      comment = get_random_name(20);
      size_t seq_size = get_random(100, 2000);
      seq = get_random_seq(seq_size);
      qual = get_random_name(seq_size);

      std::string newline = get_random(0, 1) == 1 ? "\r\n" : "\n";

      random_seqs << '@' << id << comment_spaces << comment << newline << seq
                  << newline << '+' << newline << qual << newline;

      generated_ids.push_back(id);
      generated_comments.push_back(comment);
      generated_seqs.push_back(seq);
      generated_quals.push_back(qual);
    }
    random_seqs.close();

    btllib::SeqReader random_reader(random_filename,
                                    btllib::SeqReader::Flag::LONG_MODE);
    for (i = 0; (record = random_reader.read()); i++) {
      TEST_ASSERT_EQ(record.id, generated_ids[i]);
      TEST_ASSERT_EQ(record.comment, generated_comments[i]);
      TEST_ASSERT_EQ(record.seq, generated_seqs[i]);
      TEST_ASSERT_EQ(record.qual, generated_quals[i]);
    }
    TEST_ASSERT_EQ(i, 500);

    random_reader.close();

    std::cerr << "Test larger randomly generated FASTQ file in parallel"
              << std::endl;
    std::vector<long> read_nums;

    std::vector<std::string> parallel_ids;
    std::vector<std::string> parallel_comments;
    std::vector<std::string> parallel_seqs;
    std::vector<std::string> parallel_quals;

    btllib::SeqReader random_reader2(random_filename,
                                     btllib::SeqReader::Flag::LONG_MODE);
#pragma omp parallel private(record) shared(random_reader2)
    {
      while ((record = random_reader2.read())) {
#pragma omp critical
        {
          read_nums.push_back(record.num);
          parallel_ids.push_back(record.id);
          parallel_comments.push_back(record.comment);
          parallel_seqs.push_back(record.seq);
          parallel_quals.push_back(record.qual);
        }
      }
    }

    std::sort(read_nums.begin(), read_nums.end());

    std::sort(generated_ids.begin(), generated_ids.end());
    std::sort(generated_comments.begin(), generated_comments.end());
    std::sort(generated_seqs.begin(), generated_seqs.end());
    std::sort(generated_quals.begin(), generated_quals.end());

    std::sort(parallel_ids.begin(), parallel_ids.end());
    std::sort(parallel_comments.begin(), parallel_comments.end());
    std::sort(parallel_seqs.begin(), parallel_seqs.end());
    std::sort(parallel_quals.begin(), parallel_quals.end());

    for (i = 0; i < parallel_ids.size(); i++) {
      TEST_ASSERT_EQ(read_nums[i], long(i));
      TEST_ASSERT_EQ(parallel_ids[i], generated_ids[i]);
      TEST_ASSERT_EQ(parallel_comments[i], generated_comments[i]);
      TEST_ASSERT_EQ(parallel_seqs[i], generated_seqs[i]);
      TEST_ASSERT_EQ(parallel_quals[i], generated_quals[i]);
    }
    TEST_ASSERT_EQ(i, 500);
    TEST_ASSERT_EQ(size_t(i), read_nums.size());

    std::cerr << "Test larger randomly generated FASTQ file in parallel "
                 "(different syntax)"
              << std::endl;

    read_nums.clear();

    parallel_ids.clear();
    parallel_comments.clear();
    parallel_seqs.clear();
    parallel_quals.clear();

    btllib::SeqReader random_reader3(random_filename,
                                     btllib::SeqReader::Flag::LONG_MODE);
#pragma omp parallel
    for (const auto record : random_reader3) {
#pragma omp critical
      {
        read_nums.push_back(record.num);
        parallel_ids.push_back(record.id);
        parallel_comments.push_back(record.comment);
        parallel_seqs.push_back(record.seq);
        parallel_quals.push_back(record.qual);
      }
    }

    std::sort(read_nums.begin(), read_nums.end());

    std::sort(parallel_ids.begin(), parallel_ids.end());
    std::sort(parallel_comments.begin(), parallel_comments.end());
    std::sort(parallel_seqs.begin(), parallel_seqs.end());
    std::sort(parallel_quals.begin(), parallel_quals.end());

    for (i = 0; i < parallel_ids.size(); i++) {
      TEST_ASSERT_EQ(read_nums[i], long(i));
      TEST_ASSERT_EQ(parallel_ids[i], generated_ids[i]);
      TEST_ASSERT_EQ(parallel_comments[i], generated_comments[i]);
      TEST_ASSERT_EQ(parallel_seqs[i], generated_seqs[i]);
      TEST_ASSERT_EQ(parallel_quals[i], generated_quals[i]);
    }
    TEST_ASSERT_EQ(i, 500);
    TEST_ASSERT_EQ(size_t(i), read_nums.size());

    std::remove(random_filename.c_str());
  }

  return 0;
}