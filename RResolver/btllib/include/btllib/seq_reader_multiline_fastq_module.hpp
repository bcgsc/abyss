#ifndef BTLLIB_SEQ_READER_MULTILINE_FASTQ_MODULE_HPP
#define BTLLIB_SEQ_READER_MULTILINE_FASTQ_MODULE_HPP

#include "cstring.hpp"
#include "seq.hpp"
#include <cstdlib>

namespace btllib {

/// @cond HIDDEN_SYMBOLS
class SeqReaderMultilineFastqModule
{

private:
  friend class SeqReader;

  enum class Stage
  {
    HEADER,
    SEQ,
    TRANSITION,
    SEP,
    QUAL
  };

  Stage stage = Stage::HEADER;
  CString tmp;

  static bool buffer_valid(const char* buffer, size_t size);
  template<typename ReaderType, typename RecordType>
  bool read_buffer(ReaderType& reader, RecordType& record);
  template<typename ReaderType, typename RecordType>
  bool read_transition(ReaderType& reader, RecordType& record);
  template<typename ReaderType, typename RecordType>
  bool read_file(ReaderType& reader, RecordType& record);
};

inline bool
SeqReaderMultilineFastqModule::buffer_valid(const char* buffer,
                                            const size_t size)
{
  size_t current = 0;
  unsigned char c;
  enum State
  {
    IN_HEADER_1,
    IN_HEADER_2,
    IN_SEQ,
    IN_TRANSITION,
    IN_PLUS_2,
    IN_QUAL
  };
  size_t seqlen = 0, quallen = 0;
  State state = IN_HEADER_1;
  while (current < size) {
    c = buffer[current];
    switch (state) {
      case IN_HEADER_1:
        if (c == '@') {
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
        } else if (c != '\r') {
          if (!bool(COMPLEMENTS[c])) {
            return false;
          }
          seqlen++;
        }
        break;
      case IN_TRANSITION:
        if (c == '+') {
          state = IN_PLUS_2;
          break;
        } else if (c != '\r' && !bool(COMPLEMENTS[c])) {
          return false;
        }
        seqlen++;
        state = IN_SEQ;
        break;
      case IN_PLUS_2:
        if (c == '\n') {
          state = IN_QUAL;
        }
        break;
      case IN_QUAL:
        if (quallen < seqlen) {
          if (c != '\r' && c != '\n') {
            if (c < '!' || c > '~') {
              return false;
            }
            quallen++;
          }
        } else if (c == '\n') {
          seqlen = 0;
          quallen = 0;
          state = IN_HEADER_1;
        } else if (c != '\r') {
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
SeqReaderMultilineFastqModule::read_buffer(ReaderType& reader,
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
          rtrim(record.seq);
          stage = Stage::TRANSITION;
        }
        // fall through
        case Stage::TRANSITION: {
          c = reader.getc_buffer();
          if (c == EOF) {
            return false;
          }
          reader.ungetc_buffer(c);
          if (c == '+') {
            stage = Stage::SEP;
          } else {
            stage = Stage::SEQ;
          }
          break;
        }
        case Stage::SEP: {
          if (!reader.readline_buffer_append(tmp)) {
            return false;
          }
          stage = Stage::QUAL;
          tmp.clear();
        }
        // fallthrough
        case Stage::QUAL: {
          if (!reader.readline_buffer_append(record.qual)) {
            return false;
          }
          rtrim(record.qual);
          if (record.qual.size() == record.seq.size()) {
            stage = Stage::HEADER;
            return true;
          }
          check_error(record.qual.size() > record.seq.size(),
                      "SeqReader: Multiline FASTQ reader: Quality string is "
                      "longer than sequence string.");
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
SeqReaderMultilineFastqModule::read_transition(ReaderType& reader,
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
            rtrim(record.seq);
            stage = Stage::TRANSITION;
          }
          // fall through
          case Stage::TRANSITION: {
            c = std::fgetc(reader.source);
            if (c == EOF) {
              return false;
            }
            std::ungetc(c, reader.source);
            if (c == '+') {
              stage = Stage::SEP;
            } else {
              stage = Stage::SEQ;
            }
            break;
          }
          case Stage::SEP: {
            reader.readline_file_append(tmp, reader.source);
            stage = Stage::QUAL;
            tmp.clear();
          }
          // fallthrough
          case Stage::QUAL: {
            reader.readline_file_append(record.qual, reader.source);
            rtrim(record.qual);
            if (record.qual.size() == record.seq.size()) {
              stage = Stage::HEADER;
              return true;
            }
            check_error(record.qual.size() > record.seq.size(),
                        "SeqReader: Multiline FASTQ reader: Quality string is "
                        "longer than sequence string.");
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
SeqReaderMultilineFastqModule::read_file(ReaderType& reader, RecordType& record)
{
  if (!reader.file_at_end(reader.source)) {
    reader.readline_file(record.header, reader.source);
    int c;
    reader.readline_file(record.seq, reader.source);
    rtrim(record.seq);
    for (;;) {
      c = std::fgetc(reader.source);
      check_error(c == EOF,
                  "SeqReader: Multiline FASTQ reader: Unexpected end.");
      std::ungetc(c, reader.source);
      if (c == '+') {
        reader.readline_file(tmp, reader.source);
        reader.readline_file(record.qual, reader.source);
        rtrim(record.qual);
        size_t prevlen;
        while (record.qual.size() < record.seq.size()) {
          prevlen = record.qual.size();
          reader.readline_file_append(record.qual, reader.source);
          check_error(prevlen == record.qual.size(),
                      "SeqReader: Multiline FASTQ reader: Failed to read the "
                      "quality string.");
          rtrim(record.qual);
        }
        check_error(record.qual.size() > record.seq.size(),
                    "SeqReader: Multiline FASTQ reader: Quality string is "
                    "longer than sequence string.");
        return true;
      }
      reader.readline_file_append(record.seq, reader.source);
      rtrim(record.seq);
    }
  }
  return false;
}
/// @endcond

} // namespace btllib

#endif