#ifndef BTLLIB_SEQ_WRITER_HPP
#define BTLLIB_SEQ_WRITER_HPP

#include "data_stream.hpp"
#include "seq.hpp"

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

inline SeqWriter::SeqWriter(const std::string& sink_path,
                            Format format,
                            bool append)
  : sink_path(sink_path)
  , sink(sink_path, append)
  , closed(false)
  , format(format)
  , headerchar(format == FASTA ? '>' : '@')
{}

inline void
SeqWriter::close()
{
  if (!closed) {
    sink.close();
    closed = true;
  }
}

inline void
SeqWriter::write(const std::string& id,
                 const std::string& comment,
                 const std::string& seq,
                 const std::string& qual)
{
  check_error(seq.empty(), "Attempted to write empty sequence.");
  for (const auto& c : seq) {
    if (!bool(COMPLEMENTS[(unsigned char)(c)])) {
      log_error(std::string("A sequence contains invalid IUPAC character: ") +
                c);
      std::exit(EXIT_FAILURE); // NOLINT(concurrency-mt-unsafe)
    }
  }

  std::string output;
  output.reserve(1 + id.size() + 1 + comment.size() + 1 + seq.size() + 3 +
                 qual.size() + 1);
  output += headerchar;
  if (!id.empty()) {
    output += id;
  }
  if (!comment.empty()) {
    output += " ";
    output += comment;
  }
  output += '\n';

  output += seq;
  output += '\n';

  if (format == FASTQ) {
    check_error(seq.size() != qual.size(),
                "Quality must be the same length as sequence.");
    output += "+\n";
    output += qual;
    output += '\n';
  }

  {
    std::unique_lock<std::mutex> lock(mutex);
    check_error(fwrite(output.c_str(), 1, output.size(), sink) != output.size(),
                "SeqWriter: fwrite failed.");
  }
}

} // namespace btllib

#endif