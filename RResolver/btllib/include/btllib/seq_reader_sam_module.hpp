#ifndef BTLLIB_SEQ_READER_SAM_MODULE_HPP
#define BTLLIB_SEQ_READER_SAM_MODULE_HPP

#include "btllib/cstring.hpp"
#include "btllib/process_pipeline.hpp"
#include "btllib/seq.hpp"
#include "btllib/status.hpp"

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
  loader_thread =
    std::unique_ptr<std::thread>(new std::thread([this, &reader]() {
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
          const auto ret = std::ungetc(p, reader.source);
          check_error(ret == EOF, "SeqReaderSamModule: ungetc failed.");
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
  if (!reader.file_at_end(samtools_process->out)) {
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
