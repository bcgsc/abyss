#ifndef BTLLIB_SEQ_READER_MULTILINE_FASTA_MODULE_HPP
#define BTLLIB_SEQ_READER_MULTILINE_FASTA_MODULE_HPP

#include "cstring.hpp"
#include "seq.hpp"

namespace btllib {

/// @cond HIDDEN_SYMBOLS
class SeqReaderMultilineFastaModule
{

private:
  friend class SeqReader;

  enum class Stage
  {
    HEADER,
    SEQ,
    TRANSITION,
    TRANSITION_SEQ
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
SeqReaderMultilineFastaModule::buffer_valid(const char* buffer,
                                            const size_t size)
{
  size_t current = 0;
  unsigned char c;
  enum State
  {
    IN_HEADER_1,
    IN_HEADER_2,
    IN_SEQ,
    IN_TRANSITION
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
          state = IN_TRANSITION;
        } else if (c != '\r' && !bool(COMPLEMENTS[c])) {
          return false;
        }
        break;
      case IN_TRANSITION:
        if (c == '>') {
          state = IN_HEADER_2;
          break;
        } else if (c != '\r' && !bool(COMPLEMENTS[c])) {
          return false;
        }
        state = IN_SEQ;
        break;
    }
    current++;
  }
  return true;
}

template<typename ReaderType, typename RecordType>
inline bool
SeqReaderMultilineFastaModule::read_buffer(ReaderType& reader,
                                           RecordType& record)
{
  record.header.clear();
  record.seq.clear();
  record.qual.clear();
  if (reader.buffer.start < reader.buffer.end) {
    int c;
    for (;;) {
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
          if (record.seq.back() == '\n') {
            record.seq.pop_back();
          }
          stage = Stage::TRANSITION;
        }
        // fall through
        case Stage::TRANSITION: {
          c = reader.getc_buffer();
          if (c == EOF) {
            return false;
          }
          reader.ungetc_buffer(c);
          if (c == '>') {
            stage = Stage::HEADER;
            return true;
          }
          stage = Stage::TRANSITION_SEQ;
        }
        // fall through
        case Stage::TRANSITION_SEQ: {
          if (stage == Stage::TRANSITION_SEQ) {
            if (!reader.readline_buffer_append(record.seq)) {
              return false;
            }
            if (record.seq.back() == '\n') {
              record.seq.pop_back();
            }
            stage = Stage::TRANSITION;
          }
          break;
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
SeqReaderMultilineFastaModule::read_transition(ReaderType& reader,
                                               RecordType& record)
{
  if (std::ferror(reader.source) == 0 && std::feof(reader.source) == 0) {
    const auto p = std::fgetc(reader.source);
    if (p != EOF) {
      std::ungetc(p, reader.source);
      int c;
      for (;;) {
        switch (stage) {
          case Stage::HEADER: {
            reader.readline_file_append(record.header, reader.source);
            stage = Stage::SEQ;
          }
          // fall through
          case Stage::SEQ: {
            reader.readline_file_append(record.seq, reader.source);
            if (record.seq.back() == '\n') {
              record.seq.pop_back();
            }
            stage = Stage::TRANSITION;
          }
          // fall through
          case Stage::TRANSITION: {
            c = std::fgetc(reader.source);
            if (c == EOF) {
              return false;
            }
            std::ungetc(c, reader.source);
            if (c == '>') {
              stage = Stage::HEADER;
              return true;
            }
            stage = Stage::TRANSITION_SEQ;
          }
          // fall through
          case Stage::TRANSITION_SEQ: {
            if (stage == Stage::TRANSITION_SEQ) {
              reader.readline_file_append(record.seq, reader.source);
              if (record.seq.back() == '\n') {
                record.seq.pop_back();
              }
              stage = Stage::TRANSITION;
            }
            break;
          }
          default: {
            log_error("SeqReader has entered an invalid state.");
            std::exit(EXIT_FAILURE); // NOLINT(concurrency-mt-unsafe)
          }
        }
      }
    }
  }
  return false;
}

template<typename ReaderType, typename RecordType>
inline bool
SeqReaderMultilineFastaModule::read_file(ReaderType& reader, RecordType& record)
{
  if (std::ferror(reader.source) == 0 && std::feof(reader.source) == 0) {
    reader.readline_file(record.header, reader.source);
    int c;
    reader.readline_file(record.seq, reader.source);
    if (record.seq.back() == '\n') {
      record.seq.pop_back();
    }
    for (;;) {
      c = std::fgetc(reader.source);
      if (c == EOF) {
        return true;
      }
      std::ungetc(c, reader.source);
      if (c == '>') {
        return true;
      }
      reader.readline_file_append(record.seq, reader.source);
      if (record.seq.back() == '\n') {
        record.seq.pop_back();
      }
    }
  }
  return false;
}
/// @endcond

} // namespace btllib

#endif