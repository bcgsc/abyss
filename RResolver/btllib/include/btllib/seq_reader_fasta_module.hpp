#ifndef BTLLIB_SEQ_READER_FASTA_MODULE_HPP
#define BTLLIB_SEQ_READER_FASTA_MODULE_HPP

#include "seq.hpp"

namespace btllib {

/// @cond HIDDEN_SYMBOLS
class SeqReaderFastaModule
{

private:
  friend class SeqReader;

  enum class Stage
  {
    HEADER,
    SEQ
  };

  Stage stage = Stage::HEADER;

  static bool buffer_valid(const char* buffer, size_t size);
  template<typename ReaderType, typename RecordType>
  bool read_buffer(ReaderType& reader, RecordType& record);
  template<typename ReaderType, typename RecordType>
  bool read_transition(ReaderType& reader, RecordType& record);
  template<typename ReaderType, typename RecordType>
  bool read_file(ReaderType& reader, RecordType& record);
};

inline bool
SeqReaderFastaModule::buffer_valid(const char* buffer, const size_t size)
{
  size_t current = 0;
  unsigned char c;
  enum State
  {
    IN_HEADER_1,
    IN_HEADER_2,
    IN_SEQ
  };
  State state = IN_HEADER_1;
  while (current < size) {
    c = buffer[current];
    switch (state) {
      case IN_HEADER_1:
        if (c == '>') {
          state = IN_HEADER_2;
        } else {
          return false;
        }
        break;
      case IN_HEADER_2:
        if (c == '\n') {
          state = IN_SEQ;
        }
        break;
      case IN_SEQ:
        if (c == '\n') {
          state = IN_HEADER_1;
        } else if (c != '\r' && !bool(COMPLEMENTS[c])) {
          return false;
        }
        break;
    }
    current++;
  }
  return true;
}

template<typename ReaderType, typename RecordType>
inline bool
SeqReaderFastaModule::read_buffer(ReaderType& reader, RecordType& record)
{
  record.header.clear();
  record.seq.clear();
  record.qual.clear();
  if (reader.buffer.start < reader.buffer.end) {
    switch (stage) {
      case Stage::HEADER: {
        if (!reader.readline_buffer_append(record.header)) {
          return false;
        }
        stage = Stage::SEQ;
      }
      // fall through
      case Stage::SEQ: {
        if (!reader.readline_buffer_append(record.seq)) {
          return false;
        }
        stage = Stage::HEADER;
        return true;
      }
      default: {
        log_error("SeqReader has entered an invalid state.");
        std::exit(EXIT_FAILURE); // NOLINT(concurrency-mt-unsafe)
      }
    }
  }
  return false;
}

template<typename ReaderType, typename RecordType>
inline bool
SeqReaderFastaModule::read_transition(ReaderType& reader, RecordType& record)
{
  if (std::ferror(reader.source) == 0 && std::feof(reader.source) == 0) {
    const auto p = std::fgetc(reader.source);
    if (p != EOF) {
      std::ungetc(p, reader.source);
      switch (stage) {
        case Stage::HEADER: {
          reader.readline_file_append(record.header, reader.source);
          stage = Stage::SEQ;
        }
        // fall through
        case Stage::SEQ: {
          reader.readline_file_append(record.seq, reader.source);
          stage = Stage::HEADER;
          return true;
        }
        default: {
          log_error("SeqReader has entered an invalid state.");
          std::exit(EXIT_FAILURE); // NOLINT(concurrency-mt-unsafe)
        }
      }
    }
  }
  return false;
}

template<typename ReaderType, typename RecordType>
inline bool
SeqReaderFastaModule::read_file(ReaderType& reader, RecordType& record)
{
  if (!reader.file_at_end(reader.source)) {
    reader.readline_file(record.header, reader.source);
    reader.readline_file(record.seq, reader.source);
    return true;
  }
  return false;
}
/// @endcond

} // namespace btllib

#endif