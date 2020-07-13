#ifndef BTLLIB_SEQ_READER_HPP
#define BTLLIB_SEQ_READER_HPP

#include "data_stream.hpp"
#include "index_queue.hpp"
#include "seq.hpp"
#include "status.hpp"

#include <algorithm>
#include <atomic>
#include <cassert>
#include <cctype>
#include <condition_variable>
#include <cstdio>
#include <cstring>
#include <mutex>
#include <stack>
#include <string>
#include <thread>

namespace btllib {

/** Read a FASTA, FASTQ, SAM, or GFA2 file. Threadsafe. */
class SeqReader
{
public:
  enum Flag
  {
    /** Fold lower-case characters to upper-case. */
    FOLD_CASE = 0,
    NO_FOLD_CASE = 1,
    /** Trim masked (lower case) characters from the ends of
     * sequences. */
    NO_TRIM_MASKED = 0,
    TRIM_MASKED = 2
  };

  SeqReader(const std::string& source_path, int flags = 0);
  ~SeqReader();

  void close() noexcept;

  bool flagFoldCase() const { return bool(~flags & NO_FOLD_CASE); }
  bool flagTrimMasked() const { return bool(flags & TRIM_MASKED); }

  enum Format
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

private:
  const std::string& source_path;
  DataSource source;
  unsigned flags = 0;
  Format format = UNDETERMINED; // Format of the source file
  bool closed = false;

  static const size_t DETERMINE_FORMAT_CHARS = 2048;
  static const size_t BUFFER_SIZE = DETERMINE_FORMAT_CHARS;

  char* buffer = nullptr;
  size_t buffer_start = 0;
  size_t buffer_end = 0;
  bool eof_newline_inserted = false;

  static const size_t RECORD_QUEUE_SIZE = 32;
  static const size_t RECORD_BLOCK_SIZE = 128;

  static const size_t CSTRING_DEFAULT_CAP = 4096;

  static const size_t MAX_SIMULTANEOUS_SEQREADERS = 256;

  struct CString
  {

    CString() { s[0] = '\0'; }
    CString(const CString&) = delete;
    CString(CString&& cstring) noexcept
    {
      std::swap(s, cstring.s);
      size = cstring.size;
      cstring.clear();
      std::swap(cap, cstring.cap);
    }
    CString(const std::string& str)
    {
      if (str.size() + 1 > cap) {
        cap = str.size() + 1;
        s = (char*)std::realloc((char*)s, cap); // NOLINT
      }
      size = str.size();
      memcpy(s, str.c_str(), size + 1);
    }

    CString& operator=(const CString&) = delete;
    CString& operator=(CString&& cstring) noexcept
    {
      std::swap(s, cstring.s);
      size = cstring.size;
      cstring.clear();
      std::swap(cap, cstring.cap);
      return *this;
    }
    CString& operator=(const std::string& str)
    {
      if (str.size() + 1 > cap) {
        cap = str.size() + 1;
        s = (char*)std::realloc((char*)s, cap); // NOLINT
      }
      size = str.size();
      memcpy(s, str.c_str(), size + 1);
      return *this;
    }

    ~CString() { free(s); } // NOLINT

    void clear()
    {
      s[0] = '\0';
      size = 0;
    }
    bool empty() const { return (ssize_t)size <= 0; }

    operator char*() const { return s; }

    char* s = (char*)std::malloc(CSTRING_DEFAULT_CAP); // NOLINT
    size_t size = 0;
    size_t cap = CSTRING_DEFAULT_CAP;
  };

  struct RecordCString
  {

    RecordCString() = default;
    RecordCString(const RecordCString&) = delete;
    RecordCString(RecordCString&& record) = default;

    RecordCString& operator=(const RecordCString&) = delete;
    RecordCString& operator=(RecordCString&& record) = default;

    CString header;
    CString seq;
    CString qual;
  };

  struct RecordCString2
  {

    RecordCString2() = default;
    RecordCString2(const RecordCString2&) = delete;
    RecordCString2(RecordCString2&& record) = default;

    RecordCString2& operator=(const RecordCString2&) = delete;
    RecordCString2& operator=(RecordCString2&& record) = default;

    CString header;
    std::string seq;
    CString qual;
  };

  struct RecordCString3
  {

    RecordCString3() = default;
    RecordCString3(const RecordCString3&) = delete;
    RecordCString3(RecordCString3&& record) = default;

    RecordCString3& operator=(const RecordCString3&) = delete;
    RecordCString3& operator=(RecordCString3&& record) = default;

    CString header;
    std::string seq;
    std::string qual;
  };

  CString tmp;

  std::thread* reader_thread = nullptr;
  std::thread* seq_copier_thread = nullptr;
  std::thread* qual_copier_thread = nullptr;
  std::thread* postprocessor_thread = nullptr;
  std::mutex format_mutex;
  std::condition_variable format_cv;
  std::atomic<bool> reader_end;
  RecordCString* reader_record = nullptr;
  IndexQueueSPSC<RecordCString, RECORD_QUEUE_SIZE, RECORD_BLOCK_SIZE>
    reader_queue;
  IndexQueueSPMC<RecordCString2, RECORD_QUEUE_SIZE, RECORD_BLOCK_SIZE>
    seq_copier_queue;
  IndexQueueSPMC<RecordCString3, RECORD_QUEUE_SIZE, RECORD_BLOCK_SIZE>
    qual_copier_queue;
  IndexQueueSPMC<Record, RECORD_QUEUE_SIZE, RECORD_BLOCK_SIZE>
    postprocessor_queue;

  // I am crying at this code, but until C++17 compliant compilers are
  // widespread, this cannot be a static inline variable
  static IndexQueueSPMC<Record, RECORD_QUEUE_SIZE, RECORD_BLOCK_SIZE>::Block*
  ready_records_array()
  {
    thread_local static IndexQueueSPMC<Record,
                                       RECORD_QUEUE_SIZE,
                                       RECORD_BLOCK_SIZE>::Block
      _ready_records_array[MAX_SIMULTANEOUS_SEQREADERS];
    return _ready_records_array;
  }

  // Also cry worthy
  static Record** ready_record_array()
  {
    thread_local static Record*
      _ready_record_array[MAX_SIMULTANEOUS_SEQREADERS];
    return _ready_record_array;
  }

  // Bad code bad
  static std::stack<unsigned>& recycled_ids() noexcept
  {
    static std::stack<unsigned> _recycled_ids;
    return _recycled_ids;
  }

  // ;-;
  static std::mutex& recycled_ids_mutex() noexcept
  {
    static std::mutex _recycled_ids_mutex;
    return _recycled_ids_mutex;
  };

  // :(
  static unsigned& last_id()
  {
    static unsigned _last_id = 0;
    return _last_id;
  }

  void generate_id();
  void recycle_id() const noexcept;
  unsigned id = 0;

  void determine_format();
  void start_reader();
  void start_seq_copier();
  void start_qual_copier();
  void start_postprocessor();

  bool load_buffer();

  bool is_fasta_buffer();
  bool is_fastq_buffer();
  bool is_sam_buffer();
  bool is_gfa2_buffer();

  bool readline_buffer_append(CString& s);
  void readline_file(CString& s);
  void readline_file_append(CString& s);

  int read_stage = 0;

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

  template<typename F>
  void read_from_buffer(
    F f,
    IndexQueueSPSC<RecordCString, RECORD_QUEUE_SIZE, RECORD_BLOCK_SIZE>::Block&
      records,
    size_t& counter);

  template<typename F>
  void read_transition(
    F f,
    IndexQueueSPSC<RecordCString, RECORD_QUEUE_SIZE, RECORD_BLOCK_SIZE>::Block&
      records,
    size_t& counter);

  template<typename F>
  void read_from_file(
    F f,
    IndexQueueSPSC<RecordCString, RECORD_QUEUE_SIZE, RECORD_BLOCK_SIZE>::Block&
      records,
    size_t& counter);

  void postprocess();
};

inline SeqReader::SeqReader(const std::string& source_path, int flags)
  : source_path(source_path)
  , source(source_path)
  , flags(flags)
  , reader_end(false)
{
  buffer = new char[BUFFER_SIZE];
  generate_id();
  start_seq_copier();
  start_qual_copier();
  start_postprocessor();
  {
    std::unique_lock<std::mutex> lock(format_mutex);
    start_reader();
    format_cv.wait(lock);
  }
}

inline SeqReader::~SeqReader()
{
  recycle_id();
  close();
  delete[] buffer;
  delete reader_thread;
  delete seq_copier_thread;
  delete qual_copier_thread;
  delete postprocessor_thread;
}

inline void
SeqReader::generate_id()
{
  std::unique_lock<std::mutex> lock(recycled_ids_mutex());
  if (recycled_ids().empty()) {
    id = ++last_id();
  } else {
    id = recycled_ids().top();
    recycled_ids().pop();
  }
}

inline void
SeqReader::recycle_id() const noexcept
{
  try {
    std::unique_lock<std::mutex> lock(recycled_ids_mutex());
    recycled_ids().push(id);
  } catch (const std::exception& e) {
    log_error("SeqReader id recycle error: " + std::string(e.what()));
    std::exit(EXIT_FAILURE);
  }
}

inline void
SeqReader::close() noexcept
{
  if (!closed) {
    try {
      closed = true;
      reader_end = true;
      postprocessor_queue.close();
      postprocessor_thread->join();
      qual_copier_queue.close();
      qual_copier_thread->join();
      seq_copier_queue.close();
      seq_copier_thread->join();
      reader_queue.close();
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
      fread(buffer + buffer_end, 1, BUFFER_SIZE - buffer_end, source);
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
    if (s.size >= s.cap) {
      s.cap *= 2;
      s.s = (char*)std::realloc((char*)(s.s), s.cap); // NOLINT
    }
    s.s[s.size++] = c;
  }
  if (s.size >= s.cap) {
    s.cap *= 2;
    s.s = (char*)std::realloc((char*)(s.s), s.cap); // NOLINT
  }
  s.s[s.size] = '\0';
  if (c == '\n') {
    ++buffer_start;
    return true;
  }
  return false;
}

inline void
SeqReader::readline_file(CString& s)
{
  s.size = getline(&(s.s), &(s.cap), source);
}

inline void
SeqReader::readline_file_append(CString& s)
{
  readline_file(tmp);
  if (s.size + tmp.size + 1 > s.cap) {
    s.cap = s.size + tmp.size + 1;
    s.s = (char*)std::realloc((char*)(s.s), s.cap); // NOLINT
  }
  memcpy(s.s + s.size, tmp.s, tmp.size + 1);
  s.size += tmp.size;
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
      if (tmp_string.size() + 1 > seq_reader.reader_record->header.cap) {      \
        seq_reader.reader_record->header.cap = tmp_string.size() + 1;          \
        seq_reader.reader_record->header.s =                                   \
          (char*)std::realloc((char*)(seq_reader.reader_record->header),       \
                              seq_reader.reader_record->header.cap);           \
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
      if (tmp_string.size() + 1 > seq_reader.reader_record->seq.cap) {         \
        seq_reader.reader_record->seq.cap = tmp_string.size() + 1;             \
        seq_reader.reader_record->seq.s =                                      \
          (char*)std::realloc((char*)(seq_reader.reader_record->seq.s),        \
                              seq_reader.reader_record->seq.cap);              \
      }                                                                        \
      if (tmp_string.size() + 1 > seq_reader.reader_record->qual.cap) {        \
        seq_reader.reader_record->qual.cap = tmp_string.size() + 1;            \
        seq_reader.reader_record->qual.s =                                     \
          (char*)std::realloc((char*)(seq_reader.reader_record->qual.s),       \
                              seq_reader.reader_record->qual.cap);             \
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
      if (tmp_string.size() + 1 > seq_reader.reader_record->header.cap) {      \
        seq_reader.reader_record->header.cap = tmp_string.size() + 1;          \
        seq_reader.reader_record->header.s =                                   \
          (char*)std::realloc((char*)(seq_reader.reader_record->header.s),     \
                              seq_reader.reader_record->header.cap);           \
      }                                                                        \
      seq_reader.reader_record->header = tmp_string.substr(1, pos2 - 1);       \
      for (int i = 0; i < int(SEQ) - 1; i++) {                                 \
        pos = tmp_string.find('\t', pos + 1);                                  \
      }                                                                        \
      pos2 = tmp_string.find('\t', pos + 1);                                   \
      if (pos2 == std::string::npos) {                                         \
        pos2 = tmp_string.length();                                            \
      }                                                                        \
      if (tmp_string.size() + 1 > seq_reader.reader_record->seq.cap) {         \
        seq_reader.reader_record->seq.cap = tmp_string.size() + 1;             \
        seq_reader.reader_record->seq.s =                                      \
          (char*)std::realloc((char*)(seq_reader.reader_record->seq.s),        \
                              seq_reader.reader_record->seq.cap);              \
      }                                                                        \
      seq_reader.reader_record->seq =                                          \
        tmp_string.substr(pos + 1, pos2 - pos - 1);                            \
      MIDEND_SECTION                                                           \
    }                                                                          \
    seq_reader.tmp.clear();                                                    \
    END_SECTION                                                                \
  }

struct SeqReader::read_fasta_buffer
{
  bool operator()(SeqReader& seq_reader)
  {
    switch (seq_reader.read_stage) {
      case 0: {
        if (!seq_reader.readline_buffer_append(
              seq_reader.reader_record->header)) {
          return false;
        }
        ++seq_reader.read_stage;
      }
      // fall through
      case 1: {
        if (!seq_reader.readline_buffer_append(seq_reader.reader_record->seq)) {
          return false;
        }
        seq_reader.read_stage = 0;
        return true;
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
      case 0: {
        if (!seq_reader.readline_buffer_append(
              seq_reader.reader_record->header)) {
          return false;
        }
        ++seq_reader.read_stage;
      }
      // fall through
      case 1: {
        if (!seq_reader.readline_buffer_append(seq_reader.reader_record->seq)) {
          return false;
        }
        ++seq_reader.read_stage;
      }
      // fall through
      case 2: {
        if (!seq_reader.readline_buffer_append(seq_reader.tmp)) {
          return false;
        }
        ++seq_reader.read_stage;
        seq_reader.tmp.clear();
      }
      // fall through
      case 3: {
        if (!seq_reader.readline_buffer_append(
              seq_reader.reader_record->qual)) {
          return false;
        }
        seq_reader.read_stage = 0;
        return true;
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
      case 0: {
        seq_reader.readline_file_append(seq_reader.reader_record->header);
        ++seq_reader.read_stage;
      }
      // fall through
      case 1: {
        seq_reader.readline_file_append(seq_reader.reader_record->seq);
        seq_reader.read_stage = 0;
      }
    }
  }
};

struct SeqReader::read_fastq_transition
{
  void operator()(SeqReader& seq_reader)
  {
    switch (seq_reader.read_stage) {
      case 0: {
        seq_reader.readline_file_append(seq_reader.reader_record->header);
        ++seq_reader.read_stage;
      }
      // fall through
      case 1: {
        seq_reader.readline_file_append(seq_reader.reader_record->seq);
        ++seq_reader.read_stage;
      }
      // fall through
      case 2: {
        seq_reader.readline_file_append(seq_reader.tmp);
        ++seq_reader.read_stage;
        seq_reader.tmp.clear();
      }
      // fall through
      case 3: {
        seq_reader.readline_file_append(seq_reader.reader_record->qual);
        seq_reader.read_stage = 0;
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

template<typename F>
inline void
SeqReader::read_from_buffer(
  F f,
  IndexQueueSPSC<RecordCString, RECORD_QUEUE_SIZE, RECORD_BLOCK_SIZE>::Block&
    records,
  size_t& counter)
{
  for (; buffer_start < buffer_end && !reader_end;) {
    reader_record = &(records.data[records.count]);
    if (!f(*this) || reader_record->seq.empty()) {
      break;
    }
    records.count++;
    if (records.count == RECORD_BLOCK_SIZE) {
      records.current = 0;
      records.index = counter++;
      reader_queue.write(records);
      records.current = 0;
      records.count = 0;
    }
  }
}

template<typename F>
inline void
SeqReader::read_transition(
  F f,
  IndexQueueSPSC<RecordCString, RECORD_QUEUE_SIZE, RECORD_BLOCK_SIZE>::Block&
    records,
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
        if (records.count == RECORD_BLOCK_SIZE) {
          records.current = 0;
          records.index = counter++;
          reader_queue.write(records);
          records.current = 0;
          records.count = 0;
        }
      }
    }
  }
}

template<typename F>
inline void
SeqReader::read_from_file(
  F f,
  IndexQueueSPSC<RecordCString, RECORD_QUEUE_SIZE, RECORD_BLOCK_SIZE>::Block&
    records,
  size_t& counter)
{
  for (; std::ferror(source) == 0 && std::feof(source) == 0 && !reader_end;) {
    reader_record = &(records.data[records.count]);
    f(*this);
    if (reader_record->seq.empty()) {
      break;
    }
    records.count++;
    if (records.count == RECORD_BLOCK_SIZE) {
      records.current = 0;
      records.index = counter++;
      reader_queue.write(records);
      records.current = 0;
      records.count = 0;
    }
  }
}

inline void
SeqReader::start_reader()
{
  reader_thread = new std::thread([this]() {
    {
      std::unique_lock<std::mutex> lock(format_mutex);
      determine_format();
      format_cv.notify_all();
    }

    size_t counter = 0;
    IndexQueueSPSC<RecordCString, RECORD_QUEUE_SIZE, RECORD_BLOCK_SIZE>::Block
      records;
    switch (format) {
      case FASTA: {
        read_from_buffer(read_fasta_buffer(), records, counter);
        read_transition(read_fasta_transition(), records, counter);
        read_from_file(read_fasta_file(), records, counter);
        break;
      }
      case FASTQ: {
        read_from_buffer(read_fastq_buffer(), records, counter);
        read_transition(read_fastq_transition(), records, counter);
        read_from_file(read_fastq_file(), records, counter);
        break;
      }
      case SAM: {
        read_from_buffer(read_sam_buffer(), records, counter);
        read_transition(read_sam_transition(), records, counter);
        read_from_file(read_sam_file(), records, counter);
        break;
      }
      case GFA2: {
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
    records.current = 0;
    records.index = counter++;
    size_t last_count = records.count;
    reader_queue.write(records);
    if (last_count > 0) {
      IndexQueueSPSC<RecordCString, RECORD_QUEUE_SIZE, RECORD_BLOCK_SIZE>::Block
        dummy;
      dummy.index = counter++;
      dummy.current = 0;
      dummy.count = 0;
      reader_queue.write(dummy);
    }
  });
}

inline void
SeqReader::start_seq_copier()
{
  seq_copier_thread = new std::thread([this]() {
    IndexQueueSPSC<RecordCString, RECORD_QUEUE_SIZE, RECORD_BLOCK_SIZE>::Block
      records_in;
    decltype(seq_copier_queue)::Block records_out;
    for (;;) {
      reader_queue.read(records_in);
      for (size_t i = 0; i < records_in.count; i++) {
        records_out.data[i].header = std::move(records_in.data[i].header);
        records_out.data[i].seq =
          std::string(records_in.data[i].seq.s, records_in.data[i].seq.size);
        records_out.data[i].qual = std::move(records_in.data[i].qual);
        auto& seq = records_out.data[i].seq;
        if (!seq.empty() && seq.back() == '\n') {
          seq.pop_back();
        }
      }
      records_out.count = records_in.count;
      records_out.current = records_in.current;
      records_out.index = records_in.index;
      if (records_out.count == 0) {
        seq_copier_queue.write(records_out);
        break;
      }
      seq_copier_queue.write(records_out);
    }
  });
}

inline void
SeqReader::start_qual_copier()
{
  qual_copier_thread = new std::thread([this]() {
    decltype(seq_copier_queue)::Block records_in;
    decltype(qual_copier_queue)::Block records_out;
    for (;;) {
      seq_copier_queue.read(records_in);
      for (size_t i = 0; i < records_in.count; i++) {
        records_out.data[i].header = std::move(records_in.data[i].header);
        records_out.data[i].seq = std::move(records_in.data[i].seq);
        records_out.data[i].qual =
          std::string(records_in.data[i].qual.s, records_in.data[i].qual.size);
        auto& qual = records_out.data[i].qual;
        if (!qual.empty() && qual.back() == '\n') {
          qual.pop_back();
        }
      }
      records_out.count = records_in.count;
      records_out.current = records_in.current;
      records_out.index = records_in.index;
      if (records_out.count == 0) {
        qual_copier_queue.write(records_out);
        break;
      }
      qual_copier_queue.write(records_out);
    }
  });
}

inline void
SeqReader::start_postprocessor()
{
  postprocessor_thread = new std::thread([this]() {
    decltype(qual_copier_queue)::Block records_in;
    decltype(postprocessor_queue)::Block records_out;
    for (;;) {
      qual_copier_queue.read(records_in);
      for (size_t i = 0; i < records_in.count; i++) {
        char* space = std::strstr(records_in.data[i].header, " ");
        size_t name_start = (format == FASTA || format == FASTQ) ? 1 : 0;
        if (space == nullptr) {
          records_out.data[i].name =
            std::string(records_in.data[i].header.s + name_start,
                        records_in.data[i].header.size - name_start);
          records_out.data[i].comment = "";
        } else {
          records_out.data[i].name =
            std::string(records_in.data[i].header.s + name_start,
                        space - records_in.data[i].header.s - name_start);
          records_out.data[i].comment =
            std::string(space + 1,
                        records_in.data[i].header.size -
                          (space - records_in.data[i].header.s) - 1);
        }
        records_in.data[i].header.clear();
        records_out.data[i].seq = std::move(records_in.data[i].seq);
        records_out.data[i].qual = std::move(records_in.data[i].qual);
        auto& name = records_out.data[i].name;
        auto& comment = records_out.data[i].comment;
        auto& seq = records_out.data[i].seq;
        auto& qual = records_out.data[i].qual;
        if (!name.empty() && name.back() == '\n') {
          name.pop_back();
        }
        if (!comment.empty() && comment.back() == '\n') {
          comment.pop_back();
        }
        if (flagTrimMasked()) {
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
        if (flagFoldCase()) {
          for (auto& c : seq) {
            char old = c;
            c = CAPITALS[unsigned(c)];
            if (!bool(c)) {
              log_error(std::string("A sequence contains invalid "
                                    "IUPAC character: ") +
                        old);
              std::exit(EXIT_FAILURE);
            }
          }
        }
        records_out.data[i].num = records_in.index * RECORD_BLOCK_SIZE + i;
      }
      records_out.count = records_in.count;
      records_out.current = records_in.current;
      records_out.index = records_in.index;
      if (records_out.count == 0) {
        postprocessor_queue.write(records_out);
        break;
      }
      postprocessor_queue.write(records_out);
    }
  });
}

inline SeqReader::Record
SeqReader::read()
{
  auto& ready_records = ready_records_array()[id];
  auto& ready_record = ready_record_array()[id];
  if (ready_records.count <= ready_records.current) {
    postprocessor_queue.read(ready_records);
    if (ready_records.count <= ready_records.current) {
      close();
      return Record();
    }
  }
  ready_record = &(ready_records.data[ready_records.current++]);
  return std::move(*ready_record);
}

} // namespace btllib

#endif