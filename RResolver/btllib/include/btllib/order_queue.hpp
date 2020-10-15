#ifndef BTLLIB_ORDER_QUEUE_HPP
#define BTLLIB_ORDER_QUEUE_HPP

#include <algorithm>
#include <atomic>
#include <condition_variable>
#include <mutex>
#include <string>
#include <utility>
#include <vector>

namespace btllib {

template<typename T, unsigned QUEUE_SIZE, unsigned BLOCK_SIZE>
class OrderQueue
{

public:
  struct Block
  {

    Block() {}

    Block(const Block&) = delete;

    Block(Block&& block) noexcept
      : current(block.current)
      , count(block.count)
      , num(block.num)
    {
      std::swap(data, block.data);
    }

    Block& operator=(const Block&) = delete;

    Block& operator=(Block&& block) noexcept
    {
      std::swap(data, block.data);
      current = block.current;
      count = block.count;
      num = block.num;
      return *this;
    }

    std::vector<T> data{ BLOCK_SIZE };
    size_t current = 0;
    size_t count = 0;
    size_t num = 0;
  };

  // Surrounds pieces of data in the buffer with a busy mutex
  // for exclusive access
  struct Slot
  {
    Slot() = default;
    Slot(const Slot& slot)
      : block(slot.block)
      , occupied(slot.occupied)
      , last_tenant(slot.last_tenant)
    {}

    typename OrderQueue<T, QUEUE_SIZE, BLOCK_SIZE>::Block block;
    std::mutex busy;
    bool occupied = false;
    std::condition_variable occupancy_changed;
    size_t last_tenant = -1; // Required to ensure read order
  };

  size_t elements() const { return element_count; }

  void close()
  {
    closed = true;
    for (auto& slot : this->slots) {
      slot.occupancy_changed.notify_all();
    }
  }

  bool is_closed() const { return closed; }

protected:
  std::vector<Slot> slots{ QUEUE_SIZE };
  size_t read_counter = 0;
  std::atomic<size_t> element_count{ 0 };
  std::atomic<bool> closed{ false };
};

#define ORDER_QUEUE_XPXC(SUFFIX,                                               \
                         PRE_WRITE_LOCK,                                       \
                         EXTRA_WRITE_LOCK_CONDS,                               \
                         POST_WRITE_LOCK,                                      \
                         NOTIFY_WRITE,                                         \
                         PRE_READ_LOCK,                                        \
                         EXTRA_READ_LOCK_CONDS,                                \
                         POST_READ_LOCK,                                       \
                         NOTIFY_READ,                                          \
                         MEMBERS)                                              \
  template<typename T, unsigned QUEUE_SIZE, unsigned BLOCK_SIZE>               \
  class OrderQueue##SUFFIX : public OrderQueue<T, QUEUE_SIZE, BLOCK_SIZE>      \
  {                                                                            \
                                                                               \
  public:                                                                      \
    using Block = typename OrderQueue<T, QUEUE_SIZE, BLOCK_SIZE>::Block;       \
    using Slot = typename OrderQueue<T, QUEUE_SIZE, BLOCK_SIZE>::Slot;         \
                                                                               \
    void write(Block& block)                                                   \
    {                                                                          \
      PRE_WRITE_LOCK;                                                          \
      const auto num = block.num;                                              \
      auto& target = this->slots[num % QUEUE_SIZE];                            \
      std::unique_lock<std::mutex> busy_lock(target.busy);                     \
      target.occupancy_changed.wait(busy_lock, [&] {                           \
        return (!target.occupied EXTRA_WRITE_LOCK_CONDS) || this->closed;      \
      });                                                                      \
      if (this->closed) {                                                      \
        return;                                                                \
      }                                                                        \
      POST_WRITE_LOCK; /* NOLINT */                                            \
      target.block = std::move(block);                                         \
      target.occupied = true;                                                  \
      target.occupancy_changed.NOTIFY_WRITE();                                 \
      ++(this->element_count);                                                 \
    }                                                                          \
                                                                               \
    void read(Block& block)                                                    \
    {                                                                          \
      PRE_READ_LOCK;                                                           \
      auto& target = this->slots[this->read_counter % QUEUE_SIZE];             \
      std::unique_lock<std::mutex> busy_lock(target.busy);                     \
      target.occupancy_changed.wait(busy_lock, [&] {                           \
        return (target.occupied EXTRA_READ_LOCK_CONDS) || this->closed;        \
      });                                                                      \
      if (this->closed) {                                                      \
        return;                                                                \
      }                                                                        \
      ++(this->read_counter);                                                  \
      POST_READ_LOCK;                                                          \
      block = std::move(target.block);                                         \
      target.occupied = false;                                                 \
      target.occupancy_changed.NOTIFY_READ();                                  \
      --(this->element_count);                                                 \
    }                                                                          \
                                                                               \
  private:                                                                     \
    MEMBERS; /* NOLINT */                                                      \
  };

ORDER_QUEUE_XPXC(SPSC, , , , notify_one, , , , notify_one, )
ORDER_QUEUE_XPXC(MPSC,
                 ,
                 &&(num - target.last_tenant <= QUEUE_SIZE),
                 target.last_tenant = num,
                 notify_all,
                 ,
                 ,
                 ,
                 notify_all, )
ORDER_QUEUE_XPXC(SPMC,
                 ,
                 ,
                 ,
                 notify_one,
                 std::unique_lock<std::mutex> read_lock(read_mutex),
                 ,
                 read_lock.unlock(),
                 notify_one,
                 std::mutex read_mutex)
ORDER_QUEUE_XPXC(MPMC,
                 ,
                 &&(num - target.last_tenant <= QUEUE_SIZE),
                 target.last_tenant = num,
                 notify_all,
                 std::unique_lock<std::mutex> read_lock(read_mutex),
                 ,
                 read_lock.unlock(),
                 notify_all,
                 std::mutex read_mutex)

#undef ORDER_QUEUE_XPXC

} // namespace btllib

#endif