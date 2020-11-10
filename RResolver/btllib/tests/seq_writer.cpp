#include "btllib/seq_writer.hpp"
#include "btllib/seq_reader.hpp"

#include "helpers.hpp"

#include <cstdio>
#include <fstream>

int
main()
{
  const char* names[] = { "1", "2" };
  const char* comments[] = { "comment1", "comment2" };
  const char* seqs[] = { "ACTG", "TGCA" };
  const char* quals[] = { "!@^&", "(#&$" };
  std::string random_filename;

  for (int iteration = 0; iteration < 3; iteration++) {
    std::cerr << "Iteration " << iteration + 1 << std::endl;

    // Test FASTA
    random_filename = get_random_name(64);
    std::cerr << "Test FASTA" << std::endl;
    btllib::SeqWriter writer_fasta(random_filename, btllib::SeqWriter::FASTA);
    for (int i = 0; i < 2; i++) {
      writer_fasta.write(names[i], comments[i], seqs[i], "");
    }
    writer_fasta.close();

    btllib::SeqReader reader_fasta(random_filename);
    assert(reader_fasta.get_format() == btllib::SeqReader::Format::FASTA);

    size_t i;
    btllib::SeqReader::Record record;

    i = 0;
    while ((record = reader_fasta.read())) {
      assert(record.name == names[i]);
      assert(record.comment == comments[i]);
      assert(record.seq == seqs[i]);
      assert(record.qual.empty());

      i++;
    }
    assert(i == 2);

    reader_fasta.close();
    std::remove(random_filename.c_str());

    // Test FASTQ
    random_filename = get_random_name(64) + ".bz2";
    std::cerr << "Test FASTQ" << std::endl;
    btllib::SeqWriter writer_fastq(random_filename, btllib::SeqWriter::FASTQ);
    for (int j = 0; j < 2; j++) {
      writer_fastq.write(names[j], comments[j], seqs[j], quals[j]);
    }
    writer_fastq.close();

    btllib::SeqReader reader_fastq(random_filename);
    assert(reader_fastq.get_format() == btllib::SeqReader::Format::FASTQ);

    i = 0;
    while ((record = reader_fastq.read())) {
      assert(record.name == names[i]);
      assert(record.comment == comments[i]);
      assert(record.seq == seqs[i]);
      assert(record.qual == quals[i]);

      i++;
    }
    assert(i == 2);

    reader_fastq.close();
    std::remove(random_filename.c_str());

    // Test larger randomly generated file
    std::cerr << "Test random file" << std::endl;
    std::vector<std::string> generated_names;
    std::vector<std::string> generated_comments;
    std::vector<std::string> generated_seqs;
    std::vector<std::string> generated_quals;
    random_filename = get_random_name(64) + ".gz.xz.bz2";
    btllib::SeqWriter random_seqs(
      random_filename, btllib::SeqWriter::FASTQ, false);
    for (int s = 0; s < 500; s++) {
      std::string name, comment, seq, qual;

      name = get_random_name(10);
      comment = get_random_name(20);
      size_t seq_size = get_random(100, 2000);
      seq = get_random_seq(seq_size);
      qual = get_random_name(seq_size);

      random_seqs.write(name, comment, seq, qual);

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
    std::remove(random_filename.c_str());

    std::cerr << "Test random file in parallel" << std::endl;

    auto generated_names2 = generated_names;
    auto generated_comments2 = generated_comments;
    auto generated_seqs2 = generated_seqs;
    auto generated_quals2 = generated_quals;

    auto random_filename2 = get_random_name(64) + ".gz";
    btllib::SeqWriter random_seqs2(
      random_filename2, btllib::SeqWriter::FASTQ, false);
#pragma omp parallel for
    for (int s = 0; s < 500; s++) {
      std::string name, comment, seq, qual;

#pragma omp critical
      {
        name = generated_names2.back();
        generated_names2.pop_back();
        comment = generated_comments2.back();
        generated_comments2.pop_back();
        seq = generated_seqs2.back();
        generated_seqs2.pop_back();
        qual = generated_quals2.back();
        generated_quals2.pop_back();
      }

      random_seqs2.write(name, comment, seq, qual);
    }
    random_seqs2.close();

    std::vector<std::string> parallel_names;
    std::vector<std::string> parallel_comments;
    std::vector<std::string> parallel_seqs;
    std::vector<std::string> parallel_quals;
    btllib::SeqReader random_reader2(random_filename2);
    while ((record = random_reader2.read())) {
      parallel_names.push_back(record.name);
      parallel_comments.push_back(record.comment);
      parallel_seqs.push_back(record.seq);
      parallel_quals.push_back(record.qual);
    }

    std::sort(generated_names.begin(), generated_names.end());
    std::sort(generated_comments.begin(), generated_comments.end());
    std::sort(generated_seqs.begin(), generated_seqs.end());
    std::sort(generated_quals.begin(), generated_quals.end());

    std::sort(parallel_names.begin(), parallel_names.end());
    std::sort(parallel_comments.begin(), parallel_comments.end());
    std::sort(parallel_seqs.begin(), parallel_seqs.end());
    std::sort(parallel_quals.begin(), parallel_quals.end());

    for (i = 0; i < parallel_names.size(); i++) {
      assert(parallel_names[i] == generated_names[i]);
      assert(parallel_comments[i] == generated_comments[i]);
      assert(parallel_seqs[i] == generated_seqs[i]);
      assert(parallel_quals[i] == generated_quals[i]);
    }
    assert(i == 500);

    std::remove(random_filename2.c_str());
  }

  return 0;
}