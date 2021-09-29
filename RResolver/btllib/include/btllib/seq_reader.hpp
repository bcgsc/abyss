#ifndef BTLLIB_SEQ_READER_HPP
#define BTLLIB_SEQ_READER_HPP

#include "cstring.hpp"
#include "data_stream.hpp"
#include "order_queue.hpp"
#include "seq.hpp"
#include "seq_reader_fasta_module.hpp"
#include "seq_reader_fastq_module.hpp"
#include "seq_reader_gfa2_module.hpp"
#include "seq_reader_multiline_fasta_module.hpp"
#include "seq_reader_sam_module.hpp"
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

/** Read a FASTA, FASTQ, SAM, or GFA2 file. When reading SAM files,
 * `samtools fastq` is used to convert from the SAM format to the
 * FASTQ format. Capable of reading gzip (.gz), bzip2 (.bz2), xz (.xz),
 * zip (.zip), 7zip (.7z), lrzip (.lrz), BAM (.bam) and CRAM (.cram),
 * and URL (http://, https://, ftp://) files. Threadsafe. */
class SeqReader
{
public:
  /* Has to be a struct and not an enum because:
   * 1) Non-class enums are not name qualified and can collide
   * 2) class enums can't be implicitly converted into integers
   */
  struct Flag
  {
    /** Fold all nucleotides into upper case. */
    static const unsigned FOLD_CASE = 1;
    /** Trim masked (lower case) characters from the ends of
     * sequences. */
    static const unsigned TRIM_MASKED = 2;
    /** Optimizes performance for short sequences (approx. <=5kbp) */
    static const unsigned SHORT_MODE = 4;
    /** Optimizes performance for long sequences (approx. >5kbp) */
    static const unsigned LONG_MODE = 8;
  };

  /**
   * Construct a SeqReader to read sequences from a given path.
   *
   * @param source_path Filepath to read from. Pass "-" to read from stdin.
   * @param flags Modifier flags. Specifiying either short or long mode flag is
   * mandatory; other flags are optional.
   * @param threads Maximum number of helper threads to use. Must be at least 1.
   */
  SeqReader(const std::string& source_path,
            unsigned flags,
            unsigned threads = 3);

  SeqReader(const SeqReader&) = delete;
  SeqReader(SeqReader&&) = delete;

  SeqReader& operator=(const SeqReader&) = delete;
  SeqReader& operator=(SeqReader&&) = delete;

  ~SeqReader();

  void close() noexcept;

  bool fold_case() const { return bool(flags & Flag::FOLD_CASE); }
  bool trim_masked() const { return bool(flags & Flag::TRIM_MASKED); }
  bool short_mode() const { return bool(flags & Flag::SHORT_MODE); }
  bool long_mode() const { return bool(flags & Flag::LONG_MODE); }

  enum class Format
  {
    UNDETERMINED,
    FASTA,
    MULTILINE_FASTA,
    FASTQ,
    SAM,
    GFA2,
    INVALID
  };

  friend std::ostream& operator<<(std::ostream& os, const Format f)
  {
    return os << static_cast<int32_t>(f);
  }

  Format get_format() const { return format; }

  struct Record
  {
    size_t num = -1;
    std::string id;
    std::string comment;
    std::string seq;
    std::string qual;

    operator bool() const { return !seq.empty(); }
  };

  /** Obtain next record. */
  Record read();

  static const size_t MAX_SIMULTANEOUS_SEQREADERS = 256;

  /** For range-based for loop only. */
  class RecordIterator
  {
  public:
    void operator++() { record = reader.read(); }
    bool operator!=(const RecordIterator& i)
    {
      return bool(record) || bool(i.record);
    }
    Record operator*() { return std::move(record); }
    // For wrappers
    Record next()
    {
      auto val = operator*();
      operator++();
      return val;
    }

  private:
    friend SeqReader;

    RecordIterator(SeqReader& reader, bool end)
      : reader(reader)
    {
      if (!end) {
        operator++();
      }
    }

    SeqReader& reader;
    Record record;
  };

  RecordIterator begin() { return RecordIterator(*this, false); }
  RecordIterator end() { return RecordIterator(*this, true); }

private:
  static const size_t SHORT_MODE_BUFFER_SIZE = 32;
  static const size_t SHORT_MODE_BLOCK_SIZE = 32;

  static const size_t LONG_MODE_BUFFER_SIZE = 4;
  static const size_t LONG_MODE_BLOCK_SIZE = 1;

  static const size_t FORMAT_BUFFER_SIZE = 2048;

  struct Buffer
  {

    Buffer()
      : data(FORMAT_BUFFER_SIZE)
    {}

    std::vector<char> data;
    size_t start = 0;
    size_t end = 0;
    bool eof_newline_inserted = false;
  };

  struct RecordCString
  {
    CString header;
    CString seq;
    CString qual;
  };

  const std::string& source_path;
  DataSource source;
  const unsigned flags;
  const unsigned threads;
  Format format = Format::UNDETERMINED; // Format of the source file
  std::atomic<bool> closed{ false };
  Buffer buffer;
  std::unique_ptr<std::thread> reader_thread;
  std::vector<std::unique_ptr<std::thread>> processor_threads;
  std::mutex format_mutex;
  std::condition_variable format_cv;
  std::atomic<bool> reader_end{ false };
  RecordCString* reader_record = nullptr;
  const size_t buffer_size;
  const size_t block_size;
  OrderQueueSPMC<RecordCString> cstring_queue;
  OrderQueueMPMC<Record> output_queue;
  size_t dummy_block_num = 0;
  const long id;

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

  bool load_buffer();
  void determine_format();
  void start_reader();
  void start_processors();

  CString tmp;
  bool readline_buffer_append(CString& s);
  static void readline_file(CString& s, FILE* f);
  void readline_file_append(CString& s, FILE* f);
  int getc_buffer();
  int ungetc_buffer(int c);

  inline void update_cstring_records(
    OrderQueueSPMC<RecordCString>::Block& records,
    size_t& counter);

  template<typename Module>
  void read_from_buffer(Module& module,
                        OrderQueueSPMC<RecordCString>::Block& records,
                        size_t& counter);

  template<typename Module>
  void read_transition(Module& module,
                       OrderQueueSPMC<RecordCString>::Block& records,
                       size_t& counter);

  template<typename Module>
  void read_from_file(Module& module,
                      OrderQueueSPMC<RecordCString>::Block& records,
                      size_t& counter);

  friend class SeqReaderFastaModule;
  SeqReaderFastaModule fasta_module;

  friend class SeqReaderMultilineFastaModule;
  SeqReaderMultilineFastaModule multiline_fasta_module;

  friend class SeqReaderFastqModule;
  SeqReaderFastqModule fastq_module;

  friend class SeqReaderSamModule;
  SeqReaderSamModule sam_module;

  friend class SeqReaderGfa2Module;
  SeqReaderGfa2Module gfa2_module;

  void postprocess();
};

inline SeqReader::SeqReader(const std::string& source_path,
                            const unsigned flags,
                            const unsigned threads)
  : source_path(source_path)
  , source(source_path)
  , flags(flags)
  , threads(threads)
  , buffer_size(short_mode() ? SHORT_MODE_BUFFER_SIZE : LONG_MODE_BUFFER_SIZE)
  , block_size(short_mode() ? SHORT_MODE_BLOCK_SIZE : LONG_MODE_BLOCK_SIZE)
  , cstring_queue(buffer_size, block_size)
  , output_queue(buffer_size, block_size)
  , id(++last_id())
{
  check_error(!short_mode() && !long_mode(),
              "SeqReader: no mode selected, either short or long mode flag "
              "must be provided.");
  check_error(threads == 0, "SeqReader: Number of helper threads cannot be 0.");
  start_processors();
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
  bool closed_expected = false;
  if (closed.compare_exchange_strong(closed_expected, true)) {
    try {
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
      std::exit(EXIT_FAILURE); // NOLINT(concurrency-mt-unsafe)
    }
  }
}

inline bool
SeqReader::load_buffer()
{
  buffer.start = 0;
  char last = buffer.end > 0 ? buffer.data[buffer.end - 1] : char(0);
  buffer.end = 0;
  do {
    buffer.end += fread(buffer.data.data() + buffer.end,
                        1,
                        buffer.data.size() - buffer.end,
                        source);
  } while (buffer.end < buffer.data.size() && !bool(std::feof(source)));

  if (bool(std::feof(source)) && !buffer.eof_newline_inserted) {
    if (buffer.end < buffer.data.size()) {
      if ((buffer.end == 0 && last != '\n') ||
          (buffer.end > 0 && buffer.data[buffer.end - 1] != '\n')) {
        buffer.data[buffer.end++] = '\n';
      }
      buffer.eof_newline_inserted = true;
    } else if (buffer.data[buffer.data.size() - 1] == '\n') {
      buffer.eof_newline_inserted = true;
    }
    return true;
  }
  return bool(buffer.end);
}

inline void
SeqReader::determine_format()
{
  load_buffer();
  bool empty = buffer.end - buffer.start == 1;
  check_warning(empty, std::string(source_path) + " is empty.");

  if (empty) {
    return;
  }

  auto* const buf = buffer.data.data() + buffer.start;
  const auto bufsize = buffer.end - buffer.start;

  if (fasta_module.buffer_valid(buf, bufsize)) {
    format = Format::FASTA;
  } else if (multiline_fasta_module.buffer_valid(buf, bufsize)) {
    format = Format::MULTILINE_FASTA;
  } else if (fastq_module.buffer_valid(buf, bufsize)) {
    format = Format::FASTQ;
  } else if (sam_module.buffer_valid(buf, bufsize)) {
    format = Format::SAM;
  } else if (gfa2_module.buffer_valid(buf, bufsize)) {
    format = Format::GFA2;
  } else {
    format = Format::INVALID;
    log_error(std::string(source_path) + " source file is in invalid format!");
    std::exit(EXIT_FAILURE); // NOLINT(concurrency-mt-unsafe)
  }
}

inline bool
SeqReader::readline_buffer_append(CString& s)
{
  char c = char(0);
  for (; buffer.start < buffer.end && (c = buffer.data[buffer.start]) != '\n';
       ++buffer.start) {
    if (s.s_size >= s.s_cap) {
      s.change_cap(s.s_cap * 2);
    }
    s.s[s.s_size++] = c;
  }
  if (s.s_size >= s.s_cap) {
    s.change_cap(s.s_cap * 2);
  }
  s.s[s.s_size] = '\0';
  if (c == '\n') {
    ++buffer.start;
    return true;
  }
  return false;
}

inline void
SeqReader::readline_file(CString& s, FILE* f)
{
  s.s_size = getline(&(s.s), &(s.s_cap), f);
}

inline void
SeqReader::readline_file_append(CString& s, FILE* f)
{
  readline_file(tmp, f);
  if (s.s_size + tmp.s_size + 1 > s.s_cap) {
    s.change_cap(s.s_size + tmp.s_size + 1);
  }
  memcpy(s.s + s.s_size, tmp.s, tmp.s_size + 1);
  s.s_size += tmp.s_size;
}

inline int
SeqReader::getc_buffer()
{
  if (buffer.start < buffer.end) {
    return buffer.data[buffer.start++];
  }
  return EOF;
}

inline int
SeqReader::ungetc_buffer(const int c)
{
  if (buffer.start > 0) {
    buffer.data[--buffer.start] = char(c);
    return c;
  }
  return EOF;
}

inline void
SeqReader::update_cstring_records(OrderQueueSPMC<RecordCString>::Block& records,
                                  size_t& counter)
{
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

template<typename Module>
inline void
SeqReader::read_from_buffer(Module& module,
                            OrderQueueSPMC<RecordCString>::Block& records,
                            size_t& counter)
{
  while (!reader_end) {
    reader_record = &(records.data[records.count]);
    if (!module.read_buffer(*this, *reader_record) ||
        reader_record->seq.empty()) {
      break;
    }
    update_cstring_records(records, counter);
  }
}

template<typename Module>
inline void
SeqReader::read_transition(Module& module,
                           OrderQueueSPMC<RecordCString>::Block& records,
                           size_t& counter)
{
  if (!reader_end) {
    reader_record = &(records.data[records.count]);
    module.read_transition(*this, *reader_record);
    if (!reader_record->seq.empty()) {
      update_cstring_records(records, counter);
    }
  } else if (!reader_record->seq.empty()) {
    update_cstring_records(records, counter);
  }
}

template<typename Module>
inline void
SeqReader::read_from_file(Module& module,
                          OrderQueueSPMC<RecordCString>::Block& records,
                          size_t& counter)
{
  while (!reader_end) {
    reader_record = &(records.data[records.count]);
    if (!module.read_file(*this, *reader_record) ||
        reader_record->seq.empty()) {
      break;
    }
    update_cstring_records(records, counter);
  }
}

#define BTLLIB_SEQREADER_FORMAT_READ(INPUT_FORMAT, READ_MODULE)                \
  case Format::INPUT_FORMAT: {                                                 \
    read_from_buffer(READ_MODULE, records, counter);                           \
    read_transition(READ_MODULE, records, counter);                            \
    read_from_file(READ_MODULE, records, counter);                             \
    break;                                                                     \
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
      BTLLIB_SEQREADER_FORMAT_READ(FASTA, fasta_module)
      BTLLIB_SEQREADER_FORMAT_READ(MULTILINE_FASTA, multiline_fasta_module)
      BTLLIB_SEQREADER_FORMAT_READ(FASTQ, fastq_module)
      BTLLIB_SEQREADER_FORMAT_READ(SAM, sam_module)
      BTLLIB_SEQREADER_FORMAT_READ(GFA2, gfa2_module)
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
      if (i == 0) {
        dummy_block_num = counter;
      }
      decltype(cstring_queue)::Block dummy(block_size);
      dummy.num = counter++;
      dummy.current = 0;
      dummy.count = 0;
      cstring_queue.write(dummy);
    }
  }));
}

#undef BTLLIB_SEQREADER_FORMAT_READ

inline void
SeqReader::start_processors()
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
            while (!seq.empty() && (seq.back() == '\r' || seq.back() == '\n')) {
              seq.pop_back();
            }

            records_out.data[i].qual = std::string(
              records_in.data[i].qual, records_in.data[i].qual.size());
            auto& qual = records_out.data[i].qual;
            while (!qual.empty() &&
                   (qual.back() == '\r' || qual.back() == '\n')) {
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
            size_t id_start =
              (format == Format::FASTA || format == Format::MULTILINE_FASTA ||
               format == Format::FASTQ || format == Format::SAM)
                ? 1
                : 0;

            if (first_whitespace == nullptr) {
              records_out.data[i].id =
                std::string(records_in.data[i].header + id_start,
                            records_in.data[i].header.size() - id_start);
              records_out.data[i].comment = "";
            } else {
              records_out.data[i].id = std::string(
                records_in.data[i].header + id_start,
                first_whitespace - records_in.data[i].header - id_start);
              records_out.data[i].comment = std::string(
                last_whitespace + 1,
                records_in.data[i].header.size() -
                  (last_whitespace - records_in.data[i].header) - 1);
            }
            records_in.data[i].header.clear();

            auto& id = records_out.data[i].id;
            auto& comment = records_out.data[i].comment;
            while (!id.empty() && (id.back() == '\r' || id.back() == '\n')) {
              id.pop_back();
            }
            while (!comment.empty() &&
                   (comment.back() == '\r' || comment.back() == '\n')) {
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
                  std::exit(EXIT_FAILURE); // NOLINT(concurrency-mt-unsafe)
                }
              }
            }
            records_out.data[i].num = records_in.num * block_size + i;
          }
          records_out.count = records_in.count;
          records_out.current = records_in.current;
          records_out.num = records_in.num;
          if (records_out.count == 0) {
            if (records_out.num == dummy_block_num) {
              output_queue.write(records_out);
            }
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