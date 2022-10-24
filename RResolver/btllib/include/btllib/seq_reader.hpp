#ifndef BTLLIB_SEQ_READER_HPP
#define BTLLIB_SEQ_READER_HPP

#include "btllib/cstring.hpp"
#include "btllib/data_stream.hpp"
#include "btllib/order_queue.hpp"
#include "btllib/seq.hpp"
#include "btllib/seq_reader_fasta_module.hpp"
#include "btllib/seq_reader_fastq_module.hpp"
#include "btllib/seq_reader_gfa2_module.hpp"
#include "btllib/seq_reader_multiline_fasta_module.hpp"
#include "btllib/seq_reader_multiline_fastq_module.hpp"
#include "btllib/seq_reader_sam_module.hpp"
#include "btllib/status.hpp"

#include <atomic>
#include <cctype>
#include <condition_variable>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
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
    size_t num = std::numeric_limits<size_t>::max();
    std::string id;
    std::string comment;
    std::string seq;
    std::string qual;

    operator bool() const { return !seq.empty(); }
  };

  /** Obtain next record. */
  Record read();

  /** Obtain a whole block of records. */
  OrderQueueMPMC<Record>::Block read_block();

  static const size_t MAX_SIMULTANEOUS_SEQREADERS = 256;

  /** For range-based for loop only. */
  /// @cond HIDDEN_SYMBOLS
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
  /// @endcond

  RecordIterator begin() { return RecordIterator(*this, false); }
  RecordIterator end() { return RecordIterator(*this, true); }

  size_t get_buffer_size() const { return buffer_size; }
  size_t get_block_size() const { return block_size; }

  static const size_t SHORT_MODE_BUFFER_SIZE = 32;
  static const size_t SHORT_MODE_BLOCK_SIZE = 32;

  static const size_t LONG_MODE_BUFFER_SIZE = 4;
  static const size_t LONG_MODE_BLOCK_SIZE = 1;

  static const size_t FORMAT_BUFFER_SIZE = 16384;

private:
  struct Buffer
  {

    Buffer()
      : data(FORMAT_BUFFER_SIZE)
    {
    }

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
  const std::atomic<size_t> buffer_size;
  const std::atomic<size_t> block_size;
  OrderQueueSPMC<RecordCString> cstring_queue;
  OrderQueueMPMC<Record> output_queue;
  std::atomic<size_t> dummy_block_num{ 0 };
  const long id;

  // NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
  thread_local static std::unique_ptr<decltype(output_queue)::Block>
    // NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
    ready_records_array[MAX_SIMULTANEOUS_SEQREADERS];

  // NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
  thread_local static long ready_records_owners[MAX_SIMULTANEOUS_SEQREADERS];

  // NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
  thread_local static size_t ready_records_current[MAX_SIMULTANEOUS_SEQREADERS];

  // NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
  static std::atomic<long> last_id;

  bool load_buffer();
  void determine_format();
  void start_reader();
  void start_processors();

  CString tmp;
  bool readline_buffer_append(CString& s);
  static void readline_file(CString& s, FILE* f);
  void readline_file_append(CString& s, FILE* f);
  static bool file_at_end(FILE* f);
  int getc_buffer();
  int ungetc_buffer(int c);

  void update_cstring_records(OrderQueueSPMC<RecordCString>::Block& records,
                              size_t& counter);

  /// @cond HIDDEN_SYMBOLS
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
  /// @endcond

  friend class SeqReaderFastaModule;
  SeqReaderFastaModule fasta_module;

  friend class SeqReaderMultilineFastaModule;
  SeqReaderMultilineFastaModule multiline_fasta_module;

  friend class SeqReaderFastqModule;
  SeqReaderFastqModule fastq_module;

  friend class SeqReaderMultilineFastqModule;
  SeqReaderMultilineFastqModule multiline_fastq_module;

  friend class SeqReaderSamModule;
  SeqReaderSamModule sam_module;

  friend class SeqReaderGfa2Module;
  SeqReaderGfa2Module gfa2_module;

  int module_in_use = 0;

  void postprocess();
};

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

} // namespace btllib

#endif