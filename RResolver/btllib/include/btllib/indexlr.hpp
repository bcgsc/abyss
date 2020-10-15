#ifndef BTLLIB_INDEXLR_HPP
#define BTLLIB_INDEXLR_HPP

#include "bloom_filter.hpp"
#include "nthash.hpp"
#include "order_queue.hpp"
#include "seq_reader.hpp"
#include "status.hpp"
#include "util.hpp"

#include <algorithm>
#include <atomic>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <iostream>
#include <limits>
#include <string>
#include <thread>
#include <vector>

namespace btllib {

static const unsigned BUFFER_SIZE = 256;
static const unsigned BLOCK_SIZE = 64;

// TODO: Allow multiple Indexlr objects to be instantiated (by assigning ID to
// each instance / indexing static members based on ID)
class Indexlr
{

public:
  /* Has to be a struct and not an enum because:
   * 1) Non-class enums are not name qualified and can collide
   * 2) class enums can't be implicitly converted into integers
   */
  struct Flag
  {
    static const unsigned ID = 0;
    static const unsigned NO_ID = 1;
    static const unsigned BX = 2;
    static const unsigned NO_BX = 0;
    static const unsigned SEQ = 4;
    static const unsigned NO_SEQ = 0;
    static const unsigned FILTER_IN = 8;
    static const unsigned NO_FILTER_IN = 0;
    static const unsigned FILTER_OUT = 16;
    static const unsigned NO_FILTER_OUT = 0;
  };

  bool output_id() const { return bool(~flags & Flag::NO_ID); }
  bool output_bx() const { return bool(flags & Flag::BX); }
  bool output_seq() const { return bool(flags & Flag::SEQ); }
  bool filter_in() const { return bool(flags & Flag::FILTER_IN); }
  bool filter_out() const { return bool(flags & Flag::FILTER_OUT); }

  struct Read
  {
    Read() {}

    Read(size_t num, std::string id, std::string comment, std::string seq)
      : num(num)
      , id(std::move(id))
      , comment(std::move(comment))
      , seq(std::move(seq))
    {}

    size_t num = 0;
    std::string id;
    std::string comment;
    std::string seq;
  };

  struct Minimizer
  {
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
    {}

    uint64_t min_hash, out_hash;
    size_t pos;
    bool forward;
    std::string seq;
  };

  using HashedKmer = Minimizer;

  struct Record
  {
    Record() {}

    Record(size_t num,
           std::string id,
           std::string barcode,
           std::vector<Minimizer> minimizers)
      : num(num)
      , id(std::move(id))
      , barcode(std::move(barcode))
      , minimizers(std::move(minimizers))
    {}

    size_t num = 0;
    std::string id;
    std::string barcode;
    std::vector<Minimizer> minimizers;

    operator bool() const { return !id.empty() || !barcode.empty(); }
  };

  Record get_minimizers();

  Indexlr(std::string seqfile,
          size_t k,
          size_t w,
          unsigned flags = 0,
          unsigned threads = 5,
          bool verbose = false,
          const btllib::BloomFilter& bf1 = Indexlr::dummy_bf(),
          const btllib::BloomFilter& bf2 = Indexlr::dummy_bf());

  ~Indexlr();

  static const size_t MAX_SIMULTANEOUS_INDEXLRS = 256;

private:
  static std::string extract_barcode(const std::string& id,
                                     const std::string& comment);
  std::vector<HashedKmer> hash_kmers(const std::string& seq, size_t k) const;
  static std::vector<Minimizer> minimize_hashed_kmers(
    const std::vector<HashedKmer>& hashed_kmers,
    size_t w);

  const std::string seqfile;
  const size_t k, w;
  const unsigned flags;
  const unsigned threads;
  const bool verbose;
  const unsigned id;

  static const BloomFilter& dummy_bf()
  {
    static const BloomFilter VAR;
    return VAR;
  }

  const std::reference_wrapper<const BloomFilter> bf1;
  const std::reference_wrapper<const BloomFilter> bf2;
  bool filter_in_enabled;
  bool filter_out_enabled;

  std::atomic<bool> fasta{ false };
  OrderQueueSPMC<Read, BUFFER_SIZE, BLOCK_SIZE> input_queue;
  OrderQueueMPSC<Record, BUFFER_SIZE, BLOCK_SIZE> output_queue;

  using OutputQueueType = decltype(output_queue);
  static OutputQueueType::Block* ready_blocks_array()
  {
    thread_local static decltype(
      output_queue)::Block var[MAX_SIMULTANEOUS_INDEXLRS];
    return var;
  }

  static std::atomic<unsigned>& last_id()
  {
    static std::atomic<unsigned> var(0);
    return var;
  }

  class Worker
  {
  public:
    void start() { t = std::thread(do_work, this); }
    void join() { t.join(); }

    virtual ~Worker() {}

    Worker& operator=(const Worker& worker) = delete;
    Worker& operator=(Worker&& worker) = delete;

  protected:
    Worker(Indexlr& indexlr)
      : indexlr(indexlr)
    {}

    Worker(const Worker& worker)
      : Worker(worker.indexlr)
    {}
    Worker(Worker&& worker) noexcept
      : Worker(worker.indexlr)
    {}

    Indexlr& indexlr;

    virtual void work() = 0;
    static void do_work(Worker* worker) { worker->work(); }

    std::thread t;
  };

  class InputWorker : public Worker
  {
  public:
    InputWorker(Indexlr& indexlr)
      : Worker(indexlr)
    {}

    InputWorker(const InputWorker& worker)
      : InputWorker(worker.indexlr)
    {}
    InputWorker(InputWorker&& worker) noexcept
      : InputWorker(worker.indexlr)
    {}

    InputWorker& operator=(const InputWorker& worker) = delete;
    InputWorker& operator=(InputWorker&& worker) = delete;

    void work() override;
  };

  class MinimizeWorker : public Worker
  {
  public:
    MinimizeWorker(Indexlr& indexlr)
      : Worker(indexlr)
    {}

    MinimizeWorker(const MinimizeWorker& worker)
      : MinimizeWorker(worker.indexlr)
    {}
    MinimizeWorker(MinimizeWorker&& worker) noexcept
      : MinimizeWorker(worker.indexlr)
    {}

    MinimizeWorker& operator=(const MinimizeWorker& worker) = delete;
    MinimizeWorker& operator=(MinimizeWorker&& worker) = delete;

    void work() override;
  };

  SeqReader reader;
  InputWorker input_worker;
  std::vector<MinimizeWorker> minimize_workers;
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
  , threads(threads)
  , verbose(verbose)
  , id(last_id()++)
  , bf1(bf1)
  , bf2(bf2)
  , filter_in_enabled(filter_in())
  , filter_out_enabled(filter_out())
  , reader(this->seqfile)
  , input_worker(*this)
  , minimize_workers(
      std::vector<MinimizeWorker>(threads, MinimizeWorker(*this)))
{
  input_worker.start();
  for (auto& worker : minimize_workers) {
    worker.start();
  }
}

inline Indexlr::~Indexlr()
{
  reader.close();
  for (auto& worker : minimize_workers) {
    worker.join();
  }
  input_worker.join();
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
  const static std::string BARCODE_PREFIX = "BX:Z:";
  if (starts_with(comment, BARCODE_PREFIX)) {
    const auto space_pos = comment.find(' ');
    if (space_pos != std::string::npos) {
      return comment.substr(BARCODE_PREFIX.size(),
                            space_pos - BARCODE_PREFIX.size());
    }
    return comment.substr(BARCODE_PREFIX.size());
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

inline std::vector<Indexlr::HashedKmer>
Indexlr::hash_kmers(const std::string& seq, const size_t k) const
{
  std::vector<HashedKmer> hashed_kmers;
  if (seq.size() < k) {
    return {};
  }
  hashed_kmers.reserve(seq.size() - k + 1);
  for (NtHash nh(seq, k, 2); nh.roll();) {
    hashed_kmers.emplace_back(nh.hashes()[0],
                              nh.hashes()[1],
                              nh.get_pos(),
                              nh.forward(),
                              output_seq() ? seq.substr(nh.get_pos(), k) : "");
  }
  return hashed_kmers;
}

inline std::vector<Indexlr::Minimizer>
Indexlr::minimize_hashed_kmers(
  const std::vector<Indexlr::HashedKmer>& hashed_kmers,
  const size_t w)
{
  if (hashed_kmers.size() < w) {
    return {};
  }
  std::vector<Minimizer> minimizers;
  minimizers.reserve(2 * hashed_kmers.size() / w);
  int candidate_min_pos = -1, prev_min_pos = -1;
  const auto first_it = hashed_kmers.begin();
  auto min_it = hashed_kmers.end();
  for (auto left_it = first_it; left_it < hashed_kmers.end() - w + 1;
       ++left_it) {
    const auto right_it = left_it + w;
    if (candidate_min_pos < left_it - first_it) {
      // Use of operator '<=' returns the minimum that is furthest from left.
      min_it = std::min_element(
        left_it, right_it, [](const HashedKmer& a, const HashedKmer& b) {
          return a.min_hash <= b.min_hash;
        });
    } else if (right_it[-1].min_hash <= min_it->min_hash) {
      min_it = right_it - 1;
    }
    candidate_min_pos = min_it - first_it;
    if (candidate_min_pos > prev_min_pos &&
        min_it->min_hash != std::numeric_limits<uint64_t>::max()) {
      prev_min_pos = candidate_min_pos;
      minimizers.push_back(*min_it);
    }
  }
  return minimizers;
}

inline Indexlr::Record
Indexlr::get_minimizers()
{
  auto& block = ready_blocks_array()[id % MAX_SIMULTANEOUS_INDEXLRS];
  if (block.count <= block.current) {
    output_queue.read(block);
    if (block.count <= block.current) {
      output_queue.close();
      block = decltype(output_queue)::Block();
      return Record();
    }
  }
  return std::move(block.data[block.current++]);
}

inline void
Indexlr::InputWorker::work()
{
  if (indexlr.reader.get_format() == SeqReader::Format::FASTA) {
    indexlr.fasta = true;
  } else {
    indexlr.fasta = false;
  }

  decltype(indexlr.input_queue)::Block block;
  size_t current_block_num = 0;
  SeqReader::Record record;
  Read read;
  while ((record = indexlr.reader.read())) {
    block.data[block.count++] = Read(record.num,
                                     std::move(record.name),
                                     std::move(record.comment),
                                     std::move(record.seq));
    if (block.count == BLOCK_SIZE) {
      block.num = current_block_num++;
      indexlr.input_queue.write(block);
      block.count = 0;
    }
  }
  if (block.count > 0) {
    block.num = current_block_num++;
    indexlr.input_queue.write(block);
  }
  for (unsigned i = 0; i < indexlr.threads; i++) {
    block.num = current_block_num++;
    block.current = 0;
    block.count = 0;
    indexlr.input_queue.write(block);
  }
}

inline void
Indexlr::MinimizeWorker::work()
{
  decltype(indexlr.input_queue)::Block input_block;
  decltype(indexlr.output_queue)::Block output_block;

  for (;;) {
    if (input_block.current == input_block.count) {
      if (output_block.count > 0) {
        output_block.num = input_block.num;
        indexlr.output_queue.write(output_block);
        output_block.current = 0;
        output_block.count = 0;
      }
      indexlr.input_queue.read(input_block);
    }
    if (input_block.count == 0) {
      output_block.num = input_block.num;
      output_block.current = 0;
      output_block.count = 0;
      indexlr.output_queue.write(output_block);
      break;
    }
    Read& read = input_block.data[input_block.current++];
    Record record;
    record.num = read.num;
    if (indexlr.output_id()) {
      record.id = std::move(read.id);
    }
    if (indexlr.output_bx()) {
      record.barcode = indexlr.extract_barcode(record.id, read.comment);
    }

    check_info(indexlr.verbose && read.seq.size() < indexlr.k,
               "Indexlr: skipped seq " + std::to_string(read.num) +
                 " on line " +
                 std::to_string(read.num * (indexlr.fasta ? 2 : 4) + 2) +
                 "; k (" + std::to_string(indexlr.k) + ") > seq length (" +
                 std::to_string(read.seq.size()) + ")");

    decltype(indexlr.hash_kmers(read.seq, indexlr.k)) hashed_kmers;
    if (read.seq.size() >= indexlr.k) {
      hashed_kmers = indexlr.hash_kmers(read.seq, indexlr.k);

      check_info(indexlr.verbose && indexlr.w > hashed_kmers.size(),
                 "Indexlr: skipped seq " + std::to_string(read.num) +
                   " on line " +
                   std::to_string(read.num * (indexlr.fasta ? 2 : 4) + 2) +
                   "; w (" + std::to_string(indexlr.w) + ") > # of hashes (" +
                   std::to_string(hashed_kmers.size()) + ")");
      if (indexlr.w <= hashed_kmers.size()) {
        if (indexlr.filter_in() && indexlr.filter_out()) {
          std::vector<uint64_t> tmp;
          for (auto& hk : hashed_kmers) {
            tmp = { hk.min_hash };
            if (!indexlr.bf1.get().contains(tmp) ||
                indexlr.bf2.get().contains(tmp)) {
              hk.min_hash = std::numeric_limits<uint64_t>::max();
            }
          }
        } else if (indexlr.filter_in()) {
          for (auto& hk : hashed_kmers) {
            if (!indexlr.bf1.get().contains({ hk.min_hash })) {
              hk.min_hash = std::numeric_limits<uint64_t>::max();
            }
          }
        } else if (indexlr.filter_out()) {
          for (auto& hk : hashed_kmers) {
            if (indexlr.bf1.get().contains({ hk.min_hash })) {
              hk.min_hash = std::numeric_limits<uint64_t>::max();
            }
          }
        }

        record.minimizers =
          indexlr.minimize_hashed_kmers(hashed_kmers, indexlr.w);
      } else {
        record.minimizers = {};
      }
    } else {
      record.minimizers = {};
    }

    output_block.data[output_block.count++] = std::move(record);
  }
}

} // namespace btllib

#endif