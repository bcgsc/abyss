#ifndef BTLLIB_SEQ_READER_GFA2_MODULE_HPP
#define BTLLIB_SEQ_READER_GFA2_MODULE_HPP

#include "cstring.hpp"
#include "seq.hpp"
#include <cstdlib>

namespace btllib {

/// @cond HIDDEN_SYMBOLS
class SeqReaderGfa2Module
{

private:
  friend class SeqReader;

  enum class Stage
  {
    HEADER,
    SEQ,
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
SeqReaderGfa2Module::buffer_valid(const char* buffer, const size_t size)
{
  const unsigned char specs[] = { 'H', 'S', 'F', 'E', 'G', 'O', 'U' };

  enum State
  {
    IN_ID,
    IN_ID_TAB,
    IN_REST,
    IN_IGNORED
  };

  auto is_a_spec = [&](unsigned char c) {
    bool found = false;
    for (unsigned char spec : specs) {
      if (c == spec) {
        found = true;
        break;
      }
    }
    return found;
  };

  State state = is_a_spec(buffer[0]) ? IN_ID : IN_IGNORED;
  bool has_id = false;
  size_t current = 0;
  unsigned char c;
  while (current < size) {
    c = buffer[current];
    switch (state) {
      case IN_ID:
        if (!is_a_spec(c)) {
          return false;
        }
        has_id = true;
        state = IN_ID_TAB;
        break;
      case IN_ID_TAB:
        if (c != '\t') {
          return false;
        }
        state = IN_REST;
        break;
      case IN_REST:
        if (c == '\n') {
          if (current + 1 < size) {
            state = is_a_spec(buffer[current + 1]) ? IN_ID : IN_IGNORED;
          }
        }
        break;
      case IN_IGNORED:
        if (c == '\n') {
          if (current + 1 < size) {
            state = is_a_spec(buffer[current + 1]) ? IN_ID : IN_IGNORED;
          }
        }
        break;
      default:
        break;
    }
    current++;
  }

  return has_id;
}

// NOLINTNEXTLINE
#define READ_GFA2(READLINE_SECTION, MIDEND_SECTION, END_SECTION)               \
  enum Column                                                                  \
  {                                                                            \
    S = 1,                                                                     \
    ID,                                                                        \
    LEN,                                                                       \
    SEQ                                                                        \
  };                                                                           \
  for (;;) {                                                                   \
    READLINE_SECTION                                                           \
    std::string tmp_string(tmp);                                               \
    if (tmp_string.length() > 0 && tmp_string[0] == 'S') {                     \
      size_t pos = 0, pos2 = 0;                                                \
      pos2 = tmp_string.find('\t', 1);                                         \
      if (tmp_string.size() + 1 > record.header.s_cap) {                       \
        record.header.change_cap(tmp_string.size() + 1);                       \
      }                                                                        \
      record.header = tmp_string.substr(1, pos2 - 1);                          \
      for (int i = 0; i < int(SEQ) - 1; i++) {                                 \
        pos = tmp_string.find('\t', pos + 1);                                  \
      }                                                                        \
      pos2 = tmp_string.find('\t', pos + 1);                                   \
      if (pos2 == std::string::npos) {                                         \
        pos2 = tmp_string.length();                                            \
      }                                                                        \
      if (tmp_string.size() + 1 > record.seq.s_cap) {                          \
        record.seq.change_cap(tmp_string.size() + 1);                          \
      }                                                                        \
      record.seq = tmp_string.substr(pos + 1, pos2 - pos - 1);                 \
      MIDEND_SECTION                                                           \
    }                                                                          \
    tmp.clear();                                                               \
    END_SECTION                                                                \
  }

template<typename ReaderType, typename RecordType>
inline bool
SeqReaderGfa2Module::read_buffer(ReaderType& reader, RecordType& record)
{
  log_error("GFA2 files are unsupported right now.");
  std::exit(EXIT_FAILURE); // NOLINT
  (void)reader;
  (void)record;
  return false;
  /*
  READ_GFA2(                                                           // NOLINT
    if (!reader.readline_buffer_append(                                // NOLINT
          tmp)) { return false; },                                     // NOLINT
    tmp.clear();                                                       // NOLINT
    return true;                                                       // NOLINT
    , if (reader.buffer.start >= reader.buffer.end) { return false; }) // NOLINT
  return false;
  */
}

template<typename ReaderType, typename RecordType>
inline bool
SeqReaderGfa2Module::read_transition(ReaderType& reader, RecordType& record)
{
  (void)reader;
  (void)record;
  return false;
  /*
  READ_GFA2(                                       // NOLINT
    reader.readline_file_append(tmp);              // NOLINT
    , , if (bool(feof(reader.source))) { break; }) // NOLINT
  */
}

template<typename ReaderType, typename RecordType>
inline bool
SeqReaderGfa2Module::read_file(ReaderType& reader, RecordType& record)
{
  (void)reader;
  (void)record;
  return false;
  /*
  READ_GFA2(                                       // NOLINT
    reader.readline_file(tmp);                     // NOLINT
    , , if (bool(feof(reader.source))) { break; }) // NOLINT
  */
}
/// @endcond

} // namespace btllib

#endif