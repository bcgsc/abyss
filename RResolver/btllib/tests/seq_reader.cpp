#include "../include/btllib/seq_reader.hpp"

#include "helpers.hpp"

#include <algorithm>
#include <cassert>
#include <cstdio>
#include <fstream>
#include <string>

#include <omp.h>

int
main()
{
  const char* seqs[] = { "ACTG", "TGCA" };
  const char* quals[] = { "!@^&", "(#&$" };
  std::string random_filename;

  for (int iteration = 0; iteration < 3; iteration++) {
    std::cerr << "Iteration " << iteration + 1 << std::endl;

    size_t i;
    btllib::SeqReader::Record record;

    // Test FASTA
    std::cerr << "Test FASTA" << std::endl;
    btllib::SeqReader reader_fasta("../tests/input.fa.gz.bz2.xz");
    assert(reader_fasta.get_format() == btllib::SeqReader::Format::FASTA);

    i = 0;
    while ((record = reader_fasta.read())) {
      assert(record.seq == seqs[i]);
      assert(record.qual.empty());

      i++;
    }
    assert(i == 2);

    // Test FASTQ
    std::cerr << "Test FASTQ" << std::endl;
    btllib::SeqReader reader_fastq("../tests/input.fq.tar.xz");
    assert(reader_fastq.get_format() == btllib::SeqReader::Format::FASTQ);

    i = 0;
    while ((record = reader_fastq.read())) {
      assert(record.seq == seqs[i]);
      assert(record.qual == quals[i]);

      i++;
    }
    assert(i == 2);

    // Test SAM
    std::cerr << "Test SAM" << std::endl;
    btllib::SeqReader reader_sam("../tests/input.bam");
    assert(reader_sam.get_format() == btllib::SeqReader::Format::SAM);

    i = 0;
    while ((record = reader_sam.read())) {
      assert(record.seq == seqs[i]);
      assert(record.qual == quals[i]);

      i++;
    }
    assert(i == 2);

    // Test GFA2
    std::cerr << "Test GFA2" << std::endl;
    btllib::SeqReader reader_gfa2("../tests/input.gfa2");
    assert(reader_gfa2.get_format() == btllib::SeqReader::Format::GFA2);

    i = 0;
    while ((record = reader_gfa2.read())) {
      assert(record.seq == seqs[i]);
      assert(record.qual.empty());

      i++;
    }
    assert(i == 2);

    // Test larger randomly generated file
    std::cerr << "Test random file" << std::endl;
    std::vector<std::string> generated_names;
    std::vector<std::string> generated_comments;
    std::vector<std::string> generated_seqs;
    std::vector<std::string> generated_quals;
    random_filename = get_random_name(64);
    std::ofstream random_seqs(random_filename);
    for (int s = 0; s < 500; s++) {
      std::string name, comment, seq, qual;

      name = get_random_name(10);
      comment = get_random_name(20);
      seq = get_random_sequence(200 + s);
      qual = get_random_name(200 + s);

      random_seqs << '@' << name << ' ' << comment << '\n'
                  << seq << "\n+\n"
                  << qual << '\n';

      generated_names.push_back(name);
      generated_comments.push_back(comment);
      generated_seqs.push_back(seq);
      generated_quals.push_back(qual);
    }
    random_seqs.close();

    btllib::SeqReader random_reader(random_filename);
    for (i = 0; (record = random_reader.read()); i++) {
      assert(record.name == generated_names[i]);
      assert(record.comment == generated_comments[i]);
      assert(record.seq == generated_seqs[i]);
      assert(record.qual == generated_quals[i]);
    }
    assert(i == 500);

    random_reader.close();

    std::cerr << "Test random file in parallel" << std::endl;
    std::vector<long> read_nums;

    std::vector<std::string> parallel_names;
    std::vector<std::string> parallel_comments;
    std::vector<std::string> parallel_seqs;
    std::vector<std::string> parallel_quals;

    btllib::SeqReader random_reader2(random_filename);
#pragma omp parallel private(record) shared(random_reader2)
    {
      while ((record = random_reader2.read())) {
#pragma omp critical
        {
          read_nums.push_back(record.num);
          parallel_names.push_back(record.name);
          parallel_comments.push_back(record.comment);
          parallel_seqs.push_back(record.seq);
          parallel_quals.push_back(record.qual);
        }
      }
    }

    std::sort(read_nums.begin(), read_nums.end());

    std::sort(generated_names.begin(), generated_names.end());
    std::sort(generated_comments.begin(), generated_comments.end());
    std::sort(generated_seqs.begin(), generated_seqs.end());
    std::sort(generated_quals.begin(), generated_quals.end());

    std::sort(parallel_names.begin(), parallel_names.end());
    std::sort(parallel_comments.begin(), parallel_comments.end());
    std::sort(parallel_seqs.begin(), parallel_seqs.end());
    std::sort(parallel_quals.begin(), parallel_quals.end());

    for (i = 0; i < parallel_names.size(); i++) {
      assert(read_nums[i] == long(i));
      assert(parallel_names[i] == generated_names[i]);
      assert(parallel_comments[i] == generated_comments[i]);
      assert(parallel_seqs[i] == generated_seqs[i]);
      assert(parallel_quals[i] == generated_quals[i]);
    }
    assert(i == 500);
    assert(size_t(i) == read_nums.size());

    std::remove(random_filename.c_str());
  }

  return 0;
}