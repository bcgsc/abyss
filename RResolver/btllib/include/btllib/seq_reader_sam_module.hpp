#ifndef BTLLIB_SEQ_READER_SAM_MODULE_HPP
#define BTLLIB_SEQ_READER_SAM_MODULE_HPP

#include "cstring.hpp"
#include "process_pipeline.hpp"
#include "seq.hpp"
#include "status.hpp"

#include <cstdio>
#include <memory>
#include <thread>

namespace btllib {

/// @cond HIDDEN_SYMBOLS
class SeqReaderSamModule
{

private:
  friend class SeqReader;

  enum class Stage
  {
    HEADER,
    ALIGNMENTS
  };

  std::unique_ptr<ProcessPipeline> samtools_process;
  std::unique_ptr<std::thread> loader_thread;
  CString tmp;
  static const size_t LOADER_BLOCK_SIZE = 4096;

  static bool buffer_valid(const char* buffer, size_t size);
  template<typename ReaderType, typename RecordType>
  bool read_buffer(ReaderType& reader, RecordType& record);
  template<typename ReaderType, typename RecordType>
  bool read_transition(ReaderType& reader, RecordType& record);
  template<typename ReaderType, typename RecordType>
  bool read_file(ReaderType& reader, RecordType& record);
};

inline bool
SeqReaderSamModule::buffer_valid(const char* buffer, const size_t size)
{
  enum Column
  {
    QNAME = 1,
    FLAG,
    RNAME,
    POS,
    MAPQ,
    CIGAR,
    RNEXT,
    PNEXT,
    TLEN,
    SEQ,
    QUAL
  };

  size_t current = 0;

  while (current < size && buffer[current] == '@') {
    while (current < size && buffer[current] != '\n') {
      current++;
    }
    current++;
  }

  Column column = QNAME;
  unsigned char c;
  while (current < size) {
    c = buffer[current];
    if (c == '\n') {
      break;
    }
    if (c == '\t') {
      if (current > 0 && !bool(std::isspace(buffer[current - 1]))) {
        column = Column(int(column) + 1);
      } else {
        return false;
      }
    } else {
      switch (column) {
        case QNAME:
          if (bool(std::isspace(c))) {
            return false;
          }
          break;
        case FLAG:
          if (!bool(std::isdigit(c))) {
            return false;
          }
          break;
        case RNAME:
          if (bool(std::isspace(c))) {
            return false;
          }
          break;
        case POS:
          if (!bool(std::isdigit(c))) {
            return false;
          }
          break;
        case MAPQ:
          if (!bool(std::isdigit(c))) {
            return false;
          }
          break;
        case CIGAR:
          if (bool(std::isspace(c))) {
            return false;
          }
          break;
        case RNEXT:
          if (bool(std::isspace(c))) {
            return false;
          }
          break;
        case PNEXT:
          if (!bool(std::isdigit(c))) {
            return false;
          }
          break;
        case TLEN:
          if (!bool(std::isdigit(c))) {
            return false;
          }
          break;
        case SEQ:
          if (!bool(COMPLEMENTS[c])) {
            return false;
          }
          break;
        case QUAL:
          if (bool(std::isspace(c))) {
            return false;
          }
          break;
        default:
          break;
      }
    }
    current++;
  }

  return current >= size || column >= QUAL;
}

template<typename ReaderType, typename RecordType>
inline bool
SeqReaderSamModule::read_buffer(ReaderType& reader, RecordType& record)
{
  (void)reader;
  (void)record;
  {
    ProcessPipeline version_test("samtools --version 2>/dev/stdout | head -n2");
    char* line = nullptr;
    size_t n = 0;
    std::string version = "\n";

    check_error(getline(&line, &n, version_test.out) < 2,
                "Failed to get samtools version.");
    version += line;

    line = nullptr;
    check_error(getline(&line, &n, version_test.out) < 2,
                "Failed to get samtools version.");
    version += line;
    version.pop_back();

    log_info(version);
  }
  samtools_process =
    std::unique_ptr<ProcessPipeline>(new ProcessPipeline("samtools fastq"));
  loader_thread = std::unique_ptr<std::thread>(new std::thread([&]() {
    check_error(fwrite(reader.buffer.data.data() + reader.buffer.start,
                       1,
                       reader.buffer.end - reader.buffer.start,
                       samtools_process->in) !=
                  reader.buffer.end - reader.buffer.start,
                "SeqReader SAM module: fwrite failed.");
    reader.buffer.start = reader.buffer.end;
    if (std::ferror(reader.source) == 0 && std::feof(reader.source) == 0) {
      const auto p = std::fgetc(reader.source);
      if (p != EOF) {
        std::ungetc(p, reader.source);
        while (std::ferror(reader.source) == 0 &&
               std::feof(reader.source) == 0) {
          const size_t bufsize = LOADER_BLOCK_SIZE;
          char buf[bufsize];
          size_t bytes_read = fread(buf, 1, bufsize, reader.source);
          check_error(fwrite(buf, 1, bytes_read, samtools_process->in) !=
                        bytes_read,
                      "SeqReader SAM module: fwrite failed.");
        }
      }
    }
    samtools_process->close_in();
  }));
  loader_thread->detach();
  return false;
}

template<typename ReaderType, typename RecordType>
inline bool
SeqReaderSamModule::read_transition(ReaderType& reader, RecordType& record)
{
  (void)reader;
  (void)record;
  return false;
}

template<typename ReaderType, typename RecordType>
inline bool
SeqReaderSamModule::read_file(ReaderType& reader, RecordType& record)
{
  if (std::ferror(samtools_process->out) == 0 &&
      std::feof(samtools_process->out) == 0) {
    reader.readline_file(record.header, samtools_process->out);
    reader.readline_file(record.seq, samtools_process->out);
    reader.readline_file(tmp, samtools_process->out);
    reader.readline_file(record.qual, samtools_process->out);
    return true;
  }
  return false;
}
/// @endcond

} // namespace btllib

#endif