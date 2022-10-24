#ifndef BTLLIB_SEQ_WRITER_HPP
#define BTLLIB_SEQ_WRITER_HPP

#include "btllib/data_stream.hpp"
#include "btllib/seq.hpp"

#include <cstdio>
#include <mutex>
#include <string>

namespace btllib {

/**
 * @example seq_writer.cpp
 * An example of writing a gzipped fastq file.
 */

/** Write FASTA or FASTQ sequences to a file. Capable of writing gzip (.gz),
 * bzip2 (.bz2), xz (.xz), zip (.zip), 7zip (.7z), and lrzip (.lrz) files. Add
 * the appropriate extension to the output filename to automatically compress.
 * Threadsafe. */
class SeqWriter
{

public:
  enum Format
  {
    FASTA,
    FASTQ
  };

  /**
   * Construct a SeqWriter to write sequences to a given path.
   *
   * @param source_path Filepath to write to. Pass "-" to write to stdout.
   * @param format Which format to write the output as.
   * @param append Whether to append to the target file or write anew.
   */
  SeqWriter(const std::string& sink_path,
            Format format = FASTA,
            bool append = false);

  void close();

  /**
   * Write a sequence.
   *
   * @param id Sequence ID or name.
   * @param comment Optional comment after the ID/name.
   * @param seq The sequence to write.
   * @param qual Optional quality scores, mandatory if format is FASTQ.
   */
  void write(const std::string& id,
             const std::string& comment,
             const std::string& seq,
             const std::string& qual = "");

private:
  const std::string sink_path;
  DataSink sink;
  bool closed;
  Format format;
  char headerchar;
  std::mutex mutex;
};

} // namespace btllib

#endif