#ifndef BTLLIB_SEQ_READER_HPP
#define BTLLIB_SEQ_READER_HPP

#include "cstring.hpp"
#include "data_stream.hpp"
#include "order_queue.hpp"
#include "seq.hpp"
#include "status.hpp"

#include <algorithm>
#include <atomic>
#include <cassert>
#include <cctype>
#include <condition_variable>
#include <cstdio>
#include <cstring>
#include <memory>
#include <mutex>
#include <stack>
#include <string>
#include <thread>
#include <vector>

namespace btllib {

/**
 * @example seq_reader.cpp
 * An example of reading a gzipped fastq file.
 */

/** Read a FASTA, FASTQ, SAM, or GFA2 file. Threadsafe. */
class SeqReader
{
public:
  /* Has to be a struct and not an enum because:
   * 1) Non-class enums are not name qualified and can collide
   * 2) class enums can't be implicitly converted into integers
   */
  struct Flag
  {
    /** Fold lower-case characters to upper-case. */
    static const unsigned FOLD_CASE = 0;
    static const unsigned NO_FOLD_CASE = 1;
    /** Trim masked (lower case) characters from the ends of
     * sequences. */
    static const unsigned NO_TRIM_MASKED = 0;
    static const unsigned TRIM_MASKED = 2;
  };

  SeqReader(const std::string& source_path,
            unsigned flags = 0,
            unsigned threads = 3,
            size_t buffer_size = 32,
            size_t block_size = 32);

  SeqReader(const SeqReader&) = delete;
  SeqReader(SeqReader&&) = delete;

  SeqReader& operator=(const SeqReader&) = delete;
  SeqReader& operator=(SeqReader&&) = delete;

  ~SeqReader();

  void close() noexcept;

  bool fold_case() const { return bool(~flags & Flag::NO_FOLD_CASE); }
  bool trim_masked() const { return bool(flags & Flag::TRIM_MASKED); }

  enum class Format
  {
    UNDETERMINED,
    FASTA,
    FASTQ,
    SAM,
    GFA2,
    INVALID
  };

  Format get_format() const { return format; }

  struct Record
  {
    size_t num = -1;
    std::string name;
    std::string comment;
    std::string seq;
    std::string qual;

    operator bool() const { return !seq.empty(); }
  };

  /** Read operator. */
  Record read();

  static const size_t MAX_SIMULTANEOUS_SEQREADERS = 256;

private:
  const std::string& source_path;
  DataSource source;
  const unsigned flags;
  const unsigned threads;
  Format format = Format::UNDETERMINED; // Format of the source file
  bool closed = false;

  static const size_t DETERMINE_FORMAT_CHARS = 2048;
  static const size_t BUFFER_SIZE = DETERMINE_FORMAT_CHARS;

  std::vector<char> buffer;
  size_t buffer_start = 0;
  size_t buffer_end = 0;
  bool eof_newline_inserted = false;

  struct RecordCString
  {
    CString header;
    CString seq;
    CString qual;
  };

  CString tmp;

  std::unique_ptr<std::thread> reader_thread;
  std::vector<std::unique_ptr<std::thread>> processor_threads;
  std::mutex format_mutex;
  std::condition_variable format_cv;
  std::atomic<bool> reader_end;
  RecordCString* reader_record = nullptr;
  const size_t buffer_size;
  const size_t block_size;
  OrderQueueSPMC<RecordCString> cstring_queue;
  OrderQueueMPMC<Record> output_queue;

  // I am crying at this code, but until C++17 compliant compilers are
  // widespread, this cannot be a static inline variable
  using OutputQueueType = decltype(output_queue);
  static std::unique_ptr<OutputQueueType::Block>* ready_records_array()
  {
    thread_local static std::unique_ptr<decltype(output_queue)::Block>
      var[MAX_SIMULTANEOUS_SEQREADERS];
    return var;
  }

  static long* ready_records_owners()
  {
    thread_local static long var[MAX_SIMULTANEOUS_SEQREADERS];
    return var;
  }

  // :(
  static std::atomic<long>& last_id()
  {
    static std::atomic<long> var(0);
    return var;
  }

  const long id;

  void determine_format();
  void start_reader();
  void start_processor();

  bool load_buffer();

  bool is_fasta_buffer();
  bool is_fastq_buffer();
  bool is_sam_buffer();
  bool is_gfa2_buffer();

  bool readline_buffer_append(CString& s);
  void readline_file(CString& s);
  void readline_file_append(CString& s);

  enum class ReadStage
  {
    HEADER,
    SEQ,
    SEP,
    QUAL
  };

  ReadStage read_stage = ReadStage::HEADER;

  /// @cond HIDDEN_SYMBOLS
  struct read_fasta_buffer;
  struct read_fastq_buffer;
  struct read_sam_buffer;
  struct read_gfa2_buffer;

  struct read_fasta_transition;
  struct read_fastq_transition;
  struct read_sam_transition;
  struct read_gfa2_transition;

  struct read_fasta_file;
  struct read_fastq_file;
  struct read_sam_file;
  struct read_gfa2_file;
  /// @endcond

  template<typename F>
  void read_from_buffer(F f,
                        OrderQueueSPMC<RecordCString>::Block& records,
                        size_t& counter);

  template<typename F>
  void read_transition(F f,
                       OrderQueueSPMC<RecordCString>::Block& records,
                       size_t& counter);

  template<typename F>
  void read_from_file(F f,
                      OrderQueueSPMC<RecordCString>::Block& records,
                      size_t& counter);

  void postprocess();
};

inline SeqReader::SeqReader(const std::string& source_path,
                            const unsigned flags,
                            const unsigned threads,
                            const size_t buffer_size,
                            const size_t block_size)
  : source_path(source_path)
  , source(source_path)
  , flags(flags)
  , threads(threads)
  , buffer(std::vector<char>(BUFFER_SIZE))
  , reader_end(false)
  , buffer_size(buffer_size)
  , block_size(block_size)
  , cstring_queue(buffer_size, block_size)
  , output_queue(buffer_size, block_size)
  , id(++last_id())
{
  start_processor();
  {
    std::unique_lock<std::mutex> lock(format_mutex);
    start_reader();
    format_cv.wait(lock);
  }
}

inline SeqReader::~SeqReader()
{
  close();
}

inline void
SeqReader::close() noexcept
{
  if (!closed) {
    try {
      closed = true;
      reader_end = true;
      output_queue.close();
      for (auto& pt : processor_threads) {
        pt->join();
      }
      cstring_queue.close();
      reader_thread->join();
      source.close();
    } catch (const std::system_error& e) {
      log_error("SeqReader thread join failure: " + std::string(e.what()));
      std::exit(EXIT_FAILURE);
    }
  }
}

inline bool
SeqReader::load_buffer()
{
  buffer_start = 0;
  char last = buffer_end > 0 ? buffer[buffer_end - 1] : char(0);
  buffer_end = 0;
  do {
    buffer_end +=
      fread(buffer.data() + buffer_end, 1, BUFFER_SIZE - buffer_end, source);
  } while (buffer_end < BUFFER_SIZE && !bool(std::feof(source)));

  if (bool(std::feof(source)) && !eof_newline_inserted) {
    if (buffer_end < BUFFER_SIZE) {
      if ((buffer_end == 0 && last != '\n') ||
          (buffer_end > 0 && buffer[buffer_end - 1] != '\n')) {
        buffer[buffer_end++] = '\n';
      }
      eof_newline_inserted = true;
    } else if (buffer[BUFFER_SIZE - 1] == '\n') {
      eof_newline_inserted = true;
    }
    return true;
  }
  return bool(buffer_end);
}

inline bool
SeqReader::is_fasta_buffer()
{
  size_t current = buffer_start;
  unsigned char c;
  enum State
  {
    IN_HEADER_1,
    IN_HEADER_2,
    IN_SEQ
  };
  State state = IN_HEADER_1;
  while (current < buffer_start + DETERMINE_FORMAT_CHARS &&
         current < buffer_end) {
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
        } else if (!bool(COMPLEMENTS[c])) {
          return false;
        }
        break;
    }
    current++;
  }
  return true;
}

inline bool
SeqReader::is_fastq_buffer()
{
  size_t current = buffer_start;
  unsigned char c;
  enum State
  {
    IN_HEADER_1,
    IN_HEADER_2,
    IN_SEQ,
    IN_PLUS_1,
    IN_PLUS_2,
    IN_QUAL
  };
  State state = IN_HEADER_1;
  while (current < buffer_start + DETERMINE_FORMAT_CHARS &&
         current < buffer_end) {
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
          state = IN_PLUS_1;
        } else if (!bool(COMPLEMENTS[c])) {
          return false;
        }
        break;
      case IN_PLUS_1:
        if (c == '+') {
          state = IN_PLUS_2;
        } else {
          return false;
        }
        break;
      case IN_PLUS_2:
        if (c == '\n') {
          state = IN_QUAL;
        }
        break;
      case IN_QUAL:
        if (c == '\n') {
          state = IN_HEADER_1;
        } else if (c < '!' || c > '~') {
          return false;
        }
        break;
    }
    current++;
  }
  return true;
}

inline bool
SeqReader::is_sam_buffer()
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

  size_t current = buffer_start;

  while (current < buffer_start + DETERMINE_FORMAT_CHARS &&
         current < buffer_end && buffer[current] == '@') {
    while (current < buffer_start + DETERMINE_FORMAT_CHARS &&
           current < buffer_end && buffer[current] != '\n') {
      current++;
    }
    current++;
  }

  int column = 1;
  unsigned char c;
  while (current < buffer_start + DETERMINE_FORMAT_CHARS &&
         current < buffer_end) {
    c = buffer[current];
    if (c == '\n') {
      break;
    }
    if (c == '\t') {
      if (current > 0 && !bool(std::isspace(buffer[current - 1]))) {
        column++;
      } else {
        return false;
      }
    } else {
      switch (Column(column)) {
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

  return current >= buffer_end || column >= QUAL;
}

inline bool
SeqReader::is_gfa2_buffer()
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
  size_t current = buffer_start;
  unsigned char c;
  while (current < buffer_start + DETERMINE_FORMAT_CHARS &&
         current < buffer_end) {
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
          if (current + 1 < buffer_end) {
            state = is_a_spec(buffer[current + 1]) ? IN_ID : IN_IGNORED;
          }
        }
        break;
      case IN_IGNORED:
        if (c == '\n') {
          if (current + 1 < buffer_end) {
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

inline void
SeqReader::determine_format()
{
  load_buffer();
  bool empty = buffer_end - buffer_start == 1;
  check_warning(empty, std::string(source_path) + " is empty.");

  if (empty) {
    return;
  }

  if (is_fasta_buffer()) {
    format = Format::FASTA;
  } else if (is_fastq_buffer()) {
    format = Format::FASTQ;
  } else if (is_sam_buffer()) {
    format = Format::SAM;
  } else if (is_gfa2_buffer()) {
    format = Format::GFA2;
  } else {
    format = Format::INVALID;
    log_error(std::string(source_path) + " source file is in invalid format!");
    std::exit(EXIT_FAILURE);
  }
}

inline bool
SeqReader::readline_buffer_append(CString& s)
{
  char c = char(0);
  for (; buffer_start < buffer_end && (c = buffer[buffer_start]) != '\n';
       ++buffer_start) {
    if (s.s_size >= s.s_cap) {
      s.s_cap *= 2;
      s.s = (char*)std::realloc((char*)(s.s), s.s_cap); // NOLINT
    }
    s.s[s.s_size++] = c;
  }
  if (s.s_size >= s.s_cap) {
    s.s_cap *= 2;
    s.s = (char*)std::realloc((char*)(s.s), s.s_cap); // NOLINT
  }
  s.s[s.s_size] = '\0';
  if (c == '\n') {
    ++buffer_start;
    return true;
  }
  return false;
}

inline void
SeqReader::readline_file(CString& s)
{
  s.s_size = getline(&(s.s), &(s.s_cap), source);
}

inline void
SeqReader::readline_file_append(CString& s)
{
  readline_file(tmp);
  if (s.s_size + tmp.s_size + 1 > s.s_cap) {
    s.s_cap = s.s_size + tmp.s_size + 1;
    s.s = (char*)std::realloc((char*)(s.s), s.s_cap); // NOLINT
  }
  memcpy(s.s + s.s_size, tmp.s, tmp.s_size + 1);
  s.s_size += tmp.s_size;
}

// NOLINTNEXTLINE
#define READ_SAM(READLINE_SECTION, MIDEND_SECTION, END_SECTION)                \
  enum Column                                                                  \
  {                                                                            \
    QNAME = 1,                                                                 \
    FLAG,                                                                      \
    RNAME,                                                                     \
    POS,                                                                       \
    MAPQ,                                                                      \
    CIGAR,                                                                     \
    RNEXT,                                                                     \
    PNEXT,                                                                     \
    TLEN,                                                                      \
    SEQ,                                                                       \
    QUAL                                                                       \
  };                                                                           \
  for (;;) {                                                                   \
    READLINE_SECTION                                                           \
    std::string tmp_string = seq_reader.tmp.s;                                 \
    if (tmp_string.length() > 0 && tmp_string[0] != '@') {                     \
      size_t pos = 0, pos2 = 0, pos3 = 0;                                      \
      pos2 = tmp_string.find('\t');                                            \
      if (tmp_string.size() + 1 > seq_reader.reader_record->header.s_cap) {    \
        seq_reader.reader_record->header.s_cap = tmp_string.size() + 1;        \
        seq_reader.reader_record->header.s =                                   \
          (char*)std::realloc((char*)(seq_reader.reader_record->header),       \
                              seq_reader.reader_record->header.s_cap);         \
      }                                                                        \
      seq_reader.reader_record->header = tmp_string.substr(0, pos2);           \
      for (int i = 0; i < int(SEQ) - 1; i++) {                                 \
        pos = tmp_string.find('\t', pos + 1);                                  \
      }                                                                        \
      pos2 = tmp_string.find('\t', pos + 1);                                   \
      pos3 = tmp_string.find('\t', pos2 + 1);                                  \
      if (pos3 == std::string::npos) {                                         \
        pos3 = tmp_string.length();                                            \
      }                                                                        \
      if (tmp_string.size() + 1 > seq_reader.reader_record->seq.s_cap) {       \
        seq_reader.reader_record->seq.s_cap = tmp_string.size() + 1;           \
        seq_reader.reader_record->seq.s =                                      \
          (char*)std::realloc((char*)(seq_reader.reader_record->seq.s),        \
                              seq_reader.reader_record->seq.s_cap);            \
      }                                                                        \
      if (tmp_string.size() + 1 > seq_reader.reader_record->qual.s_cap) {      \
        seq_reader.reader_record->qual.s_cap = tmp_string.size() + 1;          \
        seq_reader.reader_record->qual.s =                                     \
          (char*)std::realloc((char*)(seq_reader.reader_record->qual.s),       \
                              seq_reader.reader_record->qual.s_cap);           \
      }                                                                        \
      seq_reader.reader_record->seq =                                          \
        tmp_string.substr(pos + 1, pos2 - pos - 1);                            \
      seq_reader.reader_record->qual =                                         \
        tmp_string.substr(pos2 + 1, pos3 - pos2 - 1);                          \
      MIDEND_SECTION                                                           \
    }                                                                          \
    seq_reader.tmp.clear();                                                    \
    END_SECTION                                                                \
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
    std::string tmp_string = seq_reader.tmp.s;                                 \
    if (tmp_string.length() > 0 && tmp_string[0] == 'S') {                     \
      size_t pos = 0, pos2 = 0;                                                \
      pos2 = tmp_string.find('\t', 1);                                         \
      if (tmp_string.size() + 1 > seq_reader.reader_record->header.s_cap) {    \
        seq_reader.reader_record->header.s_cap = tmp_string.size() + 1;        \
        seq_reader.reader_record->header.s =                                   \
          (char*)std::realloc((char*)(seq_reader.reader_record->header.s),     \
                              seq_reader.reader_record->header.s_cap);         \
      }                                                                        \
      seq_reader.reader_record->header = tmp_string.substr(1, pos2 - 1);       \
      for (int i = 0; i < int(SEQ) - 1; i++) {                                 \
        pos = tmp_string.find('\t', pos + 1);                                  \
      }                                                                        \
      pos2 = tmp_string.find('\t', pos + 1);                                   \
      if (pos2 == std::string::npos) {                                         \
        pos2 = tmp_string.length();                                            \
      }                                                                        \
      if (tmp_string.size() + 1 > seq_reader.reader_record->seq.s_cap) {       \
        seq_reader.reader_record->seq.s_cap = tmp_string.size() + 1;           \
        seq_reader.reader_record->seq.s =                                      \
          (char*)std::realloc((char*)(seq_reader.reader_record->seq.s),        \
                              seq_reader.reader_record->seq.s_cap);            \
      }                                                                        \
      seq_reader.reader_record->seq =                                          \
        tmp_string.substr(pos + 1, pos2 - pos - 1);                            \
      MIDEND_SECTION                                                           \
    }                                                                          \
    seq_reader.tmp.clear();                                                    \
    END_SECTION                                                                \
  }

/// @cond HIDDEN_SYMBOLS
struct SeqReader::read_fasta_buffer
{
  bool operator()(SeqReader& seq_reader)
  {
    switch (seq_reader.read_stage) {
      case ReadStage::HEADER: {
        if (!seq_reader.readline_buffer_append(
              seq_reader.reader_record->header)) {
          return false;
        }
        seq_reader.read_stage = ReadStage::SEQ;
      }
      // fall through
      case ReadStage::SEQ: {
        if (!seq_reader.readline_buffer_append(seq_reader.reader_record->seq)) {
          return false;
        }
        seq_reader.read_stage = ReadStage::HEADER;
        return true;
      }
      default: {
        log_error("SeqReader has entered an invalid state.");
        std::exit(EXIT_FAILURE);
      }
    }
    return false;
  }
};

struct SeqReader::read_fastq_buffer
{
  bool operator()(SeqReader& seq_reader)
  {
    switch (seq_reader.read_stage) {
      case ReadStage::HEADER: {
        if (!seq_reader.readline_buffer_append(
              seq_reader.reader_record->header)) {
          return false;
        }
        seq_reader.read_stage = ReadStage::SEQ;
      }
      // fall through
      case ReadStage::SEQ: {
        if (!seq_reader.readline_buffer_append(seq_reader.reader_record->seq)) {
          return false;
        }
        seq_reader.read_stage = ReadStage::SEP;
      }
      // fall through
      case ReadStage::SEP: {
        if (!seq_reader.readline_buffer_append(seq_reader.tmp)) {
          return false;
        }
        seq_reader.read_stage = ReadStage::QUAL;
        seq_reader.tmp.clear();
      }
      // fall through
      case ReadStage::QUAL: {
        if (!seq_reader.readline_buffer_append(
              seq_reader.reader_record->qual)) {
          return false;
        }
        seq_reader.read_stage = ReadStage::HEADER;
        return true;
      }
      default: {
        log_error("SeqReader has entered an invalid state.");
        std::exit(EXIT_FAILURE);
      }
    }
    return false;
  }
};

struct SeqReader::read_sam_buffer
{
  bool operator()(SeqReader& seq_reader)
  {
    READ_SAM(                                   // NOLINT
      if (!seq_reader.readline_buffer_append(   // NOLINT
            seq_reader.tmp)) { return false; }, // NOLINT
      seq_reader.tmp.clear();                   // NOLINT
      return true;                              // NOLINT
      ,
      if (seq_reader.buffer_start >= seq_reader.buffer_end) {
        return false;
      }) // NOLINT
  }
};

struct SeqReader::read_gfa2_buffer
{
  bool operator()(SeqReader& seq_reader)
  {
    READ_GFA2(                                  // NOLINT
      if (!seq_reader.readline_buffer_append(   // NOLINT
            seq_reader.tmp)) { return false; }, // NOLINT
      seq_reader.tmp.clear();                   // NOLINT
      return true;                              // NOLINT
      ,
      if (seq_reader.buffer_start >= seq_reader.buffer_end) {
        return false;
      }) // NOLINT
  }
};

struct SeqReader::read_fasta_transition
{
  void operator()(SeqReader& seq_reader)
  {
    switch (seq_reader.read_stage) {
      case ReadStage::HEADER: {
        seq_reader.readline_file_append(seq_reader.reader_record->header);
        seq_reader.read_stage = ReadStage::SEQ;
      }
      // fall through
      case ReadStage::SEQ: {
        seq_reader.readline_file_append(seq_reader.reader_record->seq);
        seq_reader.read_stage = ReadStage::HEADER;
        return;
      }
      default: {
        log_error("SeqReader has entered an invalid state.");
        std::exit(EXIT_FAILURE);
      }
    }
  }
};

struct SeqReader::read_fastq_transition
{
  void operator()(SeqReader& seq_reader)
  {
    switch (seq_reader.read_stage) {
      case ReadStage::HEADER: {
        seq_reader.readline_file_append(seq_reader.reader_record->header);
        seq_reader.read_stage = ReadStage::SEQ;
      }
      // fall through
      case ReadStage::SEQ: {
        seq_reader.readline_file_append(seq_reader.reader_record->seq);
        seq_reader.read_stage = ReadStage::SEP;
      }
      // fall through
      case ReadStage::SEP: {
        seq_reader.readline_file_append(seq_reader.tmp);
        seq_reader.read_stage = ReadStage::QUAL;
        seq_reader.tmp.clear();
      }
      // fall through
      case ReadStage::QUAL: {
        seq_reader.readline_file_append(seq_reader.reader_record->qual);
        seq_reader.read_stage = ReadStage::HEADER;
        return;
      }
      default: {
        log_error("SeqReader has entered an invalid state.");
        std::exit(EXIT_FAILURE);
      }
    }
  }
};

struct SeqReader::read_sam_transition
{
  void operator()(SeqReader& seq_reader)
  {
    READ_SAM(                                            // NOLINT
      seq_reader.readline_file_append(seq_reader.tmp);   // NOLINT
      , , if (bool(feof(seq_reader.source))) { break; }) // NOLINT
  }
};

struct SeqReader::read_gfa2_transition
{
  void operator()(SeqReader& seq_reader)
  {
    READ_GFA2(                                           // NOLINT
      seq_reader.readline_file_append(seq_reader.tmp);   // NOLINT
      , , if (bool(feof(seq_reader.source))) { break; }) // NOLINT
  }
};

struct SeqReader::read_fasta_file
{
  void operator()(SeqReader& seq_reader)
  {
    seq_reader.readline_file(seq_reader.reader_record->header);
    seq_reader.readline_file(seq_reader.reader_record->seq);
  }
};

struct SeqReader::read_fastq_file
{
  void operator()(SeqReader& seq_reader)
  {
    seq_reader.readline_file(seq_reader.reader_record->header);
    seq_reader.readline_file(seq_reader.reader_record->seq);
    seq_reader.readline_file(seq_reader.tmp);
    seq_reader.readline_file(seq_reader.reader_record->qual);
  }
};

struct SeqReader::read_sam_file
{
  void operator()(SeqReader& seq_reader)
  {
    READ_SAM(                                            // NOLINT
      seq_reader.readline_file(seq_reader.tmp);          // NOLINT
      , , if (bool(feof(seq_reader.source))) { break; }) // NOLINT
  }
};

struct SeqReader::read_gfa2_file
{
  void operator()(SeqReader& seq_reader)
  {
    READ_GFA2(                                           // NOLINT
      seq_reader.readline_file(seq_reader.tmp);          // NOLINT
      , , if (bool(feof(seq_reader.source))) { break; }) // NOLINT
  }
};
/// @endcond

template<typename F>
inline void
SeqReader::read_from_buffer(F f,
                            OrderQueueSPMC<RecordCString>::Block& records,
                            size_t& counter)
{
  for (; buffer_start < buffer_end && !reader_end;) {
    reader_record = &(records.data[records.count]);
    reader_record->header.clear();
    reader_record->seq.clear();
    reader_record->qual.clear();
    if (!f(*this) || reader_record->seq.empty()) {
      break;
    }
    records.count++;
    if (records.count == block_size) {
      records.current = 0;
      records.num = counter++;
      cstring_queue.write(records);
      records.num = 0;
      records.current = 0;
      records.count = 0;
    }
  }
}

template<typename F>
inline void
SeqReader::read_transition(F f,
                           OrderQueueSPMC<RecordCString>::Block& records,
                           size_t& counter)
{
  if (std::ferror(source) == 0 && std::feof(source) == 0 && !reader_end) {
    int p = std::fgetc(source);
    if (p != EOF) {
      std::ungetc(p, source);
      reader_record = &(records.data[records.count]);
      f(*this);
      if (!reader_record->seq.empty()) {
        records.count++;
        if (records.count == block_size) {
          records.current = 0;
          records.num = counter++;
          cstring_queue.write(records);
          records.num = 0;
          records.current = 0;
          records.count = 0;
        }
      }
    }
  }
}

template<typename F>
inline void
SeqReader::read_from_file(F f,
                          OrderQueueSPMC<RecordCString>::Block& records,
                          size_t& counter)
{
  for (; std::ferror(source) == 0 && std::feof(source) == 0 && !reader_end;) {
    reader_record = &(records.data[records.count]);
    f(*this);
    if (reader_record->seq.empty()) {
      break;
    }
    records.count++;
    if (records.count == block_size) {
      records.current = 0;
      records.num = counter++;
      cstring_queue.write(records);
      records.num = 0;
      records.current = 0;
      records.count = 0;
    }
  }
}

inline void
SeqReader::start_reader()
{
  reader_thread = std::unique_ptr<std::thread>(new std::thread([this]() {
    {
      std::unique_lock<std::mutex> lock(format_mutex);
      determine_format();
      format_cv.notify_all();
    }

    size_t counter = 0;
    decltype(cstring_queue)::Block records(block_size);
    switch (format) {
      case Format::FASTA: {
        read_from_buffer(read_fasta_buffer(), records, counter);
        read_transition(read_fasta_transition(), records, counter);
        read_from_file(read_fasta_file(), records, counter);
        break;
      }
      case Format::FASTQ: {
        read_from_buffer(read_fastq_buffer(), records, counter);
        read_transition(read_fastq_transition(), records, counter);
        read_from_file(read_fastq_file(), records, counter);
        break;
      }
      case Format::SAM: {
        read_from_buffer(read_sam_buffer(), records, counter);
        read_transition(read_sam_transition(), records, counter);
        read_from_file(read_sam_file(), records, counter);
        break;
      }
      case Format::GFA2: {
        read_from_buffer(read_gfa2_buffer(), records, counter);
        read_transition(read_gfa2_transition(), records, counter);
        read_from_file(read_gfa2_file(), records, counter);
        break;
      }
      default: {
        break;
      }
    }

    reader_end = true;
    if (records.count > 0) {
      records.current = 0;
      records.num = counter++;
      cstring_queue.write(records);
    }
    for (unsigned i = 0; i < threads; i++) {
      decltype(cstring_queue)::Block dummy(block_size);
      dummy.num = counter++;
      dummy.current = 0;
      dummy.count = 0;
      cstring_queue.write(dummy);
    }
  }));
}

inline void
SeqReader::start_processor()
{
  processor_threads.reserve(threads);
  for (unsigned i = 0; i < threads; i++) {
    processor_threads.push_back(
      std::unique_ptr<std::thread>(new std::thread([this]() {
        decltype(cstring_queue)::Block records_in(block_size);
        decltype(output_queue)::Block records_out(block_size);
        for (;;) {
          cstring_queue.read(records_in);
          for (size_t i = 0; i < records_in.count; i++) {
            records_out.data[i].seq = std::string(
              records_in.data[i].seq, records_in.data[i].seq.size());
            auto& seq = records_out.data[i].seq;
            if (!seq.empty() && seq.back() == '\n') {
              seq.pop_back();
            }

            records_out.data[i].qual = std::string(
              records_in.data[i].qual, records_in.data[i].qual.size());
            auto& qual = records_out.data[i].qual;
            if (!qual.empty() && qual.back() == '\n') {
              qual.pop_back();
            }

            char *first_whitespace = nullptr, *last_whitespace = nullptr;
            for (size_t j = 0; j < records_in.data[i].header.size(); j++) {
              if (bool(std::isspace(records_in.data[i].header[j]))) {
                if (first_whitespace == nullptr) {
                  first_whitespace = records_in.data[i].header + j;
                }
                last_whitespace = records_in.data[i].header + j;
              } else if (last_whitespace != nullptr) {
                break;
              }
            }
            size_t name_start =
              (format == Format::FASTA || format == Format::FASTQ) ? 1 : 0;

            if (first_whitespace == nullptr) {
              records_out.data[i].name =
                std::string(records_in.data[i].header + name_start,
                            records_in.data[i].header.size() - name_start);
              records_out.data[i].comment = "";
            } else {
              records_out.data[i].name = std::string(
                records_in.data[i].header + name_start,
                first_whitespace - records_in.data[i].header - name_start);
              records_out.data[i].comment = std::string(
                last_whitespace + 1,
                records_in.data[i].header.size() -
                  (last_whitespace - records_in.data[i].header) - 1);
            }
            records_in.data[i].header.clear();

            auto& name = records_out.data[i].name;
            auto& comment = records_out.data[i].comment;
            if (!name.empty() && name.back() == '\n') {
              name.pop_back();
            }
            if (!comment.empty() && comment.back() == '\n') {
              comment.pop_back();
            }

            if (trim_masked()) {
              const auto len = seq.length();
              size_t trim_start = 0, trim_end = seq.length();
              while (trim_start <= len && bool(islower(seq[trim_start]))) {
                trim_start++;
              }
              while (trim_end > 0 && bool(islower(seq[trim_end - 1]))) {
                trim_end--;
              }
              seq.erase(trim_end);
              seq.erase(0, trim_start);
              if (!qual.empty()) {
                qual.erase(trim_end);
                qual.erase(0, trim_start);
              }
            }
            if (fold_case()) {
              for (auto& c : seq) {
                char old = c;
                c = CAPITALS[(unsigned char)(c)];
                if (!bool(c)) {
                  log_error(std::string("A sequence contains invalid "
                                        "IUPAC character: ") +
                            old);
                  std::exit(EXIT_FAILURE);
                }
              }
            }
            records_out.data[i].num = records_in.num * block_size + i;
          }
          records_out.count = records_in.count;
          records_out.current = records_in.current;
          records_out.num = records_in.num;
          if (records_out.count == 0) {
            output_queue.write(records_out);
            break;
          }
          output_queue.write(records_out);
        }
      })));
  }
}

inline SeqReader::Record
SeqReader::read()
{
  if (ready_records_owners()[id % MAX_SIMULTANEOUS_SEQREADERS] != id) {
    ready_records_array()[id % MAX_SIMULTANEOUS_SEQREADERS] =
      std::unique_ptr<decltype(output_queue)::Block>(
        new decltype(output_queue)::Block(block_size));
    ready_records_owners()[id % MAX_SIMULTANEOUS_SEQREADERS] = id;
  }
  auto& ready_records =
    *(ready_records_array()[id % MAX_SIMULTANEOUS_SEQREADERS]);
  if (ready_records.count <= ready_records.current) {
    output_queue.read(ready_records);
    if (ready_records.count <= ready_records.current) {
      close();
      ready_records = decltype(output_queue)::Block(block_size);
      return Record();
    }
  }
  return std::move(ready_records.data[ready_records.current++]);
}

} // namespace btllib

#endif