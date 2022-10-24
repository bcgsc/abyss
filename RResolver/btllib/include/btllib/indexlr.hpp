#ifndef BTLLIB_INDEXLR_HPP
#define BTLLIB_INDEXLR_HPP

#include "btllib/bloom_filter.hpp"
#include "btllib/nthash.hpp"
#include "btllib/order_queue.hpp"
#include "btllib/seq_reader.hpp"
#include "btllib/status.hpp"
#include "btllib/util.hpp"

#include <algorithm>
#include <atomic>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <thread>
#include <vector>

namespace btllib {

class Indexlr
{

public:
  /* Has to be a struct and not an enum because:
   * 1) Non-class enums are not name qualified and can collide
   * 2) class enums can't be implicitly converted into integers
   */
  struct Flag
  {
    /** Do not include read ID information. May improve performance. */
    static const unsigned NO_ID = 1;
    /** Include barcode along with minimizer information. */
    static const unsigned BX = 2;
    /** Include read sequence along with minimizer information. */
    static const unsigned SEQ = 4;
    /** Only include minimizers found in the first Bloom filter argument.
     */
    static const unsigned FILTER_IN = 8;
    /** Filter out minimizers found in a Bloom filter argument. If FILTER_IN is
     * NOT enabled, then use the first Bloom filter argument. If both FILTER_IN
     * and FILTER_OUT flags are enabled, FILTER_IN uses the first Bloom filter
     * argument and FILTER_OUT uses the second. */
    static const unsigned FILTER_OUT = 16;
    /** Optimizes performance for short sequences (approx. <=5kbp) */
    static const unsigned SHORT_MODE = 32;
    /** Optimizes performance for long sequences (approx. >5kbp) */
    static const unsigned LONG_MODE = 64;
  };

  bool output_id() const { return bool(~flags & Flag::NO_ID); }
  bool output_bx() const { return bool(flags & Flag::BX); }
  bool output_seq() const { return bool(flags & Flag::SEQ); }
  bool filter_in() const { return bool(flags & Flag::FILTER_IN); }
  bool filter_out() const { return bool(flags & Flag::FILTER_OUT); }
  bool short_mode() const { return bool(flags & Flag::SHORT_MODE); }
  bool long_mode() const { return bool(flags & Flag::LONG_MODE); }

  struct Minimizer
  {
    Minimizer() = default;

    Minimizer(uint64_t min_hash,
              uint64_t out_hash,
              size_t pos,
              bool forward,
              std::string seq)
      : min_hash(min_hash)
      , out_hash(out_hash)
      , pos(pos)
      , forward(forward)
      , seq(std::move(seq))
    {
    }

    uint64_t min_hash = 0, out_hash = 0;
    size_t pos = 0;
    bool forward = false;
    std::string seq;
  };

  using HashedKmer = Minimizer;

  struct Record
  {
    Record() {}

    Record(size_t num,
           std::string id,
           std::string barcode,
           size_t readlen,
           std::vector<Minimizer> minimizers)
      : num(num)
      , id(std::move(id))
      , barcode(std::move(barcode))
      , readlen(readlen)
      , minimizers(std::move(minimizers))
    {
    }

    size_t num = 0;
    std::string id;
    std::string barcode;
    size_t readlen = 0;
    std::vector<Minimizer> minimizers;

    operator bool() const
    {
      return !id.empty() || !barcode.empty() || !minimizers.empty();
    }
  };

  /**
   * Read the next Indexlr record, containing read
   * information and a list of minimizers.
   */
  Record read();

  /**
   * Construct Indexlr to calculate minimizers from sequences at the given path.
   *
   * @param seqfile Filepath to read sequences from. Pass "-" to read from
   * stdin.
   * @param k k-mer size for the minimizer.
   * @param w window size when selecting minimizers.
   * @param flags Modifier flags. Specifiying either short or long mode flag is
   * mandatory; other flags are optional.
   * @param threads Maximum number of processing threads to use. Must be at
   * least 1.
   * @param verbose Whether to output informational messages during processing.
   * @param bf1 A Bloom filter to use for either filtering minimizers in or out,
   * if one of the flags is enabled. If both are enabled, bf1 is used for
   * filtering in.
   * @param bf2 A second Bloom filter to use when both filtering minimizers in
   * and out is enabled. bf2 is used for filtering out.
   */
  Indexlr(std::string seqfile,
          size_t k,
          size_t w,
          unsigned flags = 0,
          unsigned threads = 5,
          bool verbose = false,
          const btllib::BloomFilter& bf1 = Indexlr::dummy_bf(),
          const btllib::BloomFilter& bf2 = Indexlr::dummy_bf());

  ~Indexlr();

  void close() noexcept;

  static const size_t MAX_SIMULTANEOUS_INDEXLRS = 256;

  /** For range-based for loop only. */
  /// @cond HIDDEN_SYMBOLS
  class RecordIterator
  {
  public:
    void operator++() { record = indexlr.read(); }
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
    friend Indexlr;

    RecordIterator(Indexlr& indexlr, bool end)
      : indexlr(indexlr)
    {
      if (!end) {
        operator++();
      }
    }

    Indexlr& indexlr;
    Record record;
  };
  /// @endcond

  RecordIterator begin() { return RecordIterator(*this, false); }
  RecordIterator end() { return RecordIterator(*this, true); }

private:
  static std::string extract_barcode(const std::string& id,
                                     const std::string& comment);
  static void filter_hashed_kmer(Indexlr::HashedKmer& hk,
                                 bool filter_in,
                                 bool filter_out,
                                 const BloomFilter& filter_in_bf,
                                 const BloomFilter& filter_out_bf);
  static void calc_minimizer(
    const std::vector<Indexlr::HashedKmer>& hashed_kmers_buffer,
    const Indexlr::Minimizer*& min_current,
    size_t idx,
    ssize_t& min_idx_left,
    ssize_t& min_idx_right,
    ssize_t& min_pos_prev,
    size_t w,
    std::vector<Indexlr::Minimizer>& minimizers);
  std::vector<Minimizer> minimize(const std::string& seq) const;

  const std::string seqfile;
  const size_t k, w;
  const unsigned flags;
  const bool verbose;
  const long id;
  std::atomic<bool> closed{ false };

  static const BloomFilter& dummy_bf()
  {
    static const BloomFilter var;
    return var;
  }

  const std::reference_wrapper<const BloomFilter> filter_in_bf;
  const std::reference_wrapper<const BloomFilter> filter_out_bf;
  bool filter_in_enabled;
  bool filter_out_enabled;

  SeqReader reader;
  OrderQueueMPSC<Record> output_queue;

  using OutputQueueType = decltype(output_queue);
  static std::unique_ptr<OutputQueueType::Block>* ready_blocks_array()
  {
    thread_local static std::unique_ptr<decltype(output_queue)::Block>
      var[MAX_SIMULTANEOUS_INDEXLRS];
    return var;
  }

  static long* ready_blocks_owners()
  {
    thread_local static long var[MAX_SIMULTANEOUS_INDEXLRS] = { 0 };
    return var;
  }

  static size_t* ready_blocks_current()
  {
    thread_local static size_t var[MAX_SIMULTANEOUS_INDEXLRS] = { 0 };
    return var;
  }

  static std::atomic<long>& last_id()
  {
    static std::atomic<long> var(0);
    return var;
  }

  class Worker
  {
  public:
    void start() { t = std::thread(do_work, this); }
    void join() { t.join(); }
    void set_id(const int id) { this->id = id; }

    Worker& operator=(const Worker& worker) = delete;
    Worker& operator=(Worker&& worker) = delete;

    Worker(Indexlr& indexlr)
      : indexlr(indexlr)
    {
    }
    Worker(const Worker& worker)
      : Worker(worker.indexlr)
    {
    }
    Worker(Worker&& worker) noexcept
      : Worker(worker.indexlr)
    {
    }

  private:
    void work();
    static void do_work(Worker* worker) { worker->work(); }

    int id = -1;
    Indexlr& indexlr;
    std::thread t;
  };

  std::vector<Worker> workers;
  Barrier end_barrier;
  std::mutex last_block_num_mutex;
  uint64_t last_block_num = 0;
  bool last_block_num_valid = false;
};

inline Indexlr::Indexlr(std::string seqfile,
                        const size_t k,
                        const size_t w,
                        const unsigned flags,
                        const unsigned threads,
                        const bool verbose,
                        const BloomFilter& bf1,
                        const BloomFilter& bf2)
  : seqfile(std::move(seqfile))
  , k(k)
  , w(w)
  , flags(flags)
  , verbose(verbose)
  , id(++last_id())
  , filter_in_bf(filter_in() ? bf1 : Indexlr::dummy_bf())
  , filter_out_bf(filter_out() ? filter_in() ? bf2 : bf1 : Indexlr::dummy_bf())
  , filter_in_enabled(filter_in())
  , filter_out_enabled(filter_out())
  , reader(this->seqfile,
           short_mode() ? SeqReader::Flag::SHORT_MODE
                        : SeqReader::Flag::LONG_MODE)
  , output_queue(reader.get_buffer_size(), reader.get_block_size())
  , workers(std::vector<Worker>(threads, Worker(*this)))
  , end_barrier(threads)
{
  check_error(!short_mode() && !long_mode(),
              "Indexlr: no mode selected, either short or long mode flag must "
              "be provided.");
  check_error(short_mode() && long_mode(),
              "Indexlr: short and long mode are mutually exclusive.");
  check_error(threads == 0,
              "Indexlr: Number of processing threads cannot be 0.");
  int id_counter = 0;
  for (auto& worker : workers) {
    worker.set_id(id_counter++);
    worker.start();
  }
}

inline Indexlr::~Indexlr()
{
  close();
}

inline void
Indexlr::close() noexcept
{
  bool closed_expected = false;
  if (closed.compare_exchange_strong(closed_expected, true)) {
    try {
      reader.close();
      output_queue.close();
      for (auto& worker : workers) {
        worker.join();
      }
    } catch (const std::system_error& e) {
      log_error("Indexlr thread join failure: " + std::string(e.what()));
      std::exit(EXIT_FAILURE); // NOLINT(concurrency-mt-unsafe)
    }
  }
}

// Minimerize a sequence: Find the minimizers of a vector of hash values
// representing a sequence.
/* Algorithm
v is a vector of non-negative integers
w is the window size
Invariants
    0 <  w <= v.size() - 1
    0 <= l <= r <= v.size() - 1
Initial conditions
    M    = NIL       Final set of minimizers, empty initially
    min  = -1        Minimum element
    i    = -1        Index of minimum element
    prev = -1        Index of previous minimum element
    l    = 0         Index of left end of window
    r    = l + w - 1 Index of right end of window
Computation
At each window, if the previous minimum is out of scope, find the new,
right-most, minimum or else, check with only the right-most element to determine
if that is the new minimum. A minimizer is added to the final vector only if
it's index has changed. for each window of v bounded by [l, r] if (i < l) i =
index of minimum element in [l, r], furthest from l. else if (v[r] <= v[i]) i =
r min = v[i] if (i != prev) { prev = i M <- M + m
    }
    l = l + 1        Move window's left bound by one element
    r = l + w - 1    Set window's right bound
}*/

inline std::string
Indexlr::extract_barcode(const std::string& id, const std::string& comment)
{
  const static std::string barcode_prefix = "BX:Z:";
  if (startswith(comment, barcode_prefix)) {
    const auto space_pos = comment.find(' ');
    if (space_pos != std::string::npos) {
      return comment.substr(barcode_prefix.size(),
                            space_pos - barcode_prefix.size());
    }
    return comment.substr(barcode_prefix.size());
  }
  const auto pound_pos = id.find('#');
  if (pound_pos != std::string::npos) {
    const auto slash_pos = id.find('/');
    if (slash_pos > pound_pos) {
      return id.substr(pound_pos + 1, slash_pos - (pound_pos + 1));
    }
  }
  return "NA";
}

inline void
Indexlr::filter_hashed_kmer(Indexlr::HashedKmer& hk,
                            bool filter_in,
                            bool filter_out,
                            const BloomFilter& filter_in_bf,
                            const BloomFilter& filter_out_bf)
{
  if (filter_in && filter_out) {
    std::vector<uint64_t> tmp;
    tmp = { hk.min_hash };
    if (!filter_in_bf.contains(tmp) || filter_out_bf.contains(tmp)) {
      hk.min_hash = std::numeric_limits<uint64_t>::max();
    }
  } else if (filter_in) {
    if (!filter_in_bf.contains({ hk.min_hash })) {
      hk.min_hash = std::numeric_limits<uint64_t>::max();
    }
  } else if (filter_out) {
    if (filter_out_bf.contains({ hk.min_hash })) {
      hk.min_hash = std::numeric_limits<uint64_t>::max();
    }
  }
}

inline void
Indexlr::calc_minimizer(
  const std::vector<Indexlr::HashedKmer>& hashed_kmers_buffer,
  const Indexlr::Minimizer*& min_current,
  const size_t idx,
  ssize_t& min_idx_left,
  ssize_t& min_idx_right,
  ssize_t& min_pos_prev,
  const size_t w,
  std::vector<Indexlr::Minimizer>& minimizers)
{
  min_idx_left = ssize_t(idx + 1 - w);
  min_idx_right = ssize_t(idx + 1);
  const auto& min_left =
    hashed_kmers_buffer[min_idx_left % hashed_kmers_buffer.size()];
  const auto& min_right =
    hashed_kmers_buffer[(min_idx_right - 1) % hashed_kmers_buffer.size()];

  if (min_current == nullptr || min_current->pos < min_left.pos) {
    min_current = &min_left;
    // Use of operator '<=' returns the minimum that is furthest from left.
    for (ssize_t i = min_idx_left; i < min_idx_right; i++) {
      const auto& min_i = hashed_kmers_buffer[i % hashed_kmers_buffer.size()];
      if (min_i.min_hash <= min_current->min_hash) {
        min_current = &min_i;
      }
    }
  } else if (min_right.min_hash <= min_current->min_hash) {
    min_current = &min_right;
  }
  if (ssize_t(min_current->pos) > min_pos_prev &&
      min_current->min_hash != std::numeric_limits<uint64_t>::max()) {
    min_pos_prev = ssize_t(min_current->pos);
    minimizers.push_back(*min_current);
  }
}

inline std::vector<Indexlr::Minimizer>
Indexlr::minimize(const std::string& seq) const
{
  if ((k > seq.size()) || (w > seq.size() - k + 1)) {
    return {};
  }
  std::vector<Minimizer> minimizers;
  minimizers.reserve(2 * (seq.size() - k + 1) / w);
  std::vector<HashedKmer> hashed_kmers_buffer(w + 1);
  ssize_t min_idx_left, min_idx_right, min_pos_prev = -1;
  const Minimizer* min_current = nullptr;
  size_t idx = 0;
  for (NtHash nh(seq, 2, k); nh.roll(); ++idx) {
    auto& hk = hashed_kmers_buffer[idx % hashed_kmers_buffer.size()];

    hk = HashedKmer(nh.hashes()[0],
                    nh.hashes()[1],
                    nh.get_pos(),
                    nh.forward(),
                    output_seq() ? seq.substr(nh.get_pos(), k) : "");

    filter_hashed_kmer(
      hk, filter_in(), filter_out(), filter_in_bf.get(), filter_out_bf.get());

    if (idx + 1 >= w) {
      calc_minimizer(hashed_kmers_buffer,
                     min_current,
                     idx,
                     min_idx_left,
                     min_idx_right,
                     min_pos_prev,
                     w,
                     minimizers);
    }
  }
  return minimizers;
}

inline Indexlr::Record
Indexlr::read()
{
  if (ready_blocks_owners()[id % MAX_SIMULTANEOUS_INDEXLRS] != id) {
    ready_blocks_array()[id % MAX_SIMULTANEOUS_INDEXLRS] =
      std::unique_ptr<decltype(output_queue)::Block>(
        new decltype(output_queue)::Block(reader.get_block_size()));
    ready_blocks_owners()[id % MAX_SIMULTANEOUS_INDEXLRS] = id;
    ready_blocks_current()[id % MAX_SIMULTANEOUS_INDEXLRS] = 0;
  }
  auto& block = *(ready_blocks_array()[id % MAX_SIMULTANEOUS_INDEXLRS]);
  auto& current = ready_blocks_current()[id % MAX_SIMULTANEOUS_INDEXLRS];
  if (current >= block.count) {
    block.count = 0;
    output_queue.read(block);
    if (block.count == 0) {
      output_queue.close();
      block = decltype(output_queue)::Block(reader.get_block_size());
      return Record();
    }
    current = 0;
  }
  return std::move(block.data[current++]);
}

inline void
Indexlr::Worker::work()
{
  decltype(indexlr.output_queue)::Block output_block(
    indexlr.reader.get_block_size());
  uint64_t last_block_num = 0;
  bool last_block_num_valid = false;
  for (;;) {
    auto input_block = indexlr.reader.read_block();
    if (input_block.count == 0) {
      break;
    }

    output_block.num = input_block.num;
    for (size_t idx = 0; idx < input_block.count; idx++) {
      Record record;
      auto& reader_record = input_block.data[idx];
      record.num = reader_record.num;
      if (indexlr.output_id()) {
        record.id = std::move(reader_record.id);
      }
      if (indexlr.output_bx()) {
        record.barcode =
          indexlr.extract_barcode(record.id, reader_record.comment);
      }
      record.readlen = reader_record.seq.size();

      check_info(indexlr.verbose && indexlr.k > record.readlen,
                 "Indexlr: skipped seq " + std::to_string(record.num) +
                   " on line " +
                   std::to_string(record.num * (indexlr.reader.get_format() ==
                                                    SeqReader::Format::FASTA
                                                  ? 2
                                                  : 4) +
                                  2) +
                   "; k (" + std::to_string(indexlr.k) + ") > seq length (" +
                   std::to_string(record.readlen) + ")");

      check_info(indexlr.verbose && indexlr.w > record.readlen - indexlr.k + 1,
                 "Indexlr: skipped seq " + std::to_string(record.num) +
                   " on line " +
                   std::to_string(record.num * (indexlr.reader.get_format() ==
                                                    SeqReader::Format::FASTA
                                                  ? 2
                                                  : 4) +
                                  2) +
                   "; w (" + std::to_string(indexlr.w) + ") > # of hashes (" +
                   std::to_string(record.readlen - indexlr.k + 1) + ")");

      if (indexlr.k <= record.readlen &&
          indexlr.w <= record.readlen - indexlr.k + 1) {
        record.minimizers = indexlr.minimize(reader_record.seq);
      } else {
        record.minimizers = {};
      }

      output_block.data[output_block.count++] = std::move(record);
    }
    if (output_block.count > 0) {
      last_block_num = output_block.num;
      last_block_num_valid = true;
      indexlr.output_queue.write(output_block);
      output_block.count = 0;
    }
  }
  if (last_block_num_valid) {
    std::unique_lock<std::mutex> lock(indexlr.last_block_num_mutex);
    indexlr.last_block_num = std::max(indexlr.last_block_num, last_block_num);
    indexlr.last_block_num_valid = true;
    lock.unlock();
  }
  indexlr.end_barrier.wait();
  if (last_block_num_valid && indexlr.last_block_num_valid &&
      last_block_num == indexlr.last_block_num) {
    output_block.num = last_block_num + 1;
    indexlr.output_queue.write(output_block);
  } else if (!indexlr.last_block_num_valid && id == 0) {
    output_block.num = 0;
    indexlr.output_queue.write(output_block);
  }
}

} // namespace btllib

#endif