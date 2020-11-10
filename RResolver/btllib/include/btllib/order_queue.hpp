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

template<typename T>
class OrderQueue
{

public:
  struct Block
  {

    Block(const size_t block_size)
      : data(block_size)
    {}

    Block(const Block& block) = default;

    Block(Block&& block) noexcept
      : current(block.current)
      , count(block.count)
      , num(block.num)
    {
      std::swap(data, block.data);
      block.current = 0;
      block.count = 0;
      block.num = 0;
    }

    Block& operator=(const Block& block) = default;

    Block& operator=(Block&& block) noexcept
    {
      std::swap(data, block.data);
      current = block.current;
      count = block.count;
      num = block.num;
      block.current = 0;
      block.count = 0;
      block.num = 0;
      return *this;
    }

    std::vector<T> data;
    size_t current = 0;
    size_t count = 0;
    size_t num = 0;
  };

  // Surrounds pieces of data in the buffer with a busy mutex
  // for exclusive access
  struct Slot
  {
    Slot(size_t block_size)
      : block(block_size)
    {}
    Slot(const Slot& slot)
      : block(slot.block)
      , occupied(slot.occupied)
      , last_tenant(slot.last_tenant)
    {}
    Slot(Slot&& slot) noexcept
      : block(slot.block)
      , occupied(slot.occupied)
      , last_tenant(slot.last_tenant)
    {}

    Slot& operator=(const Slot& slot)
    {
      if (this == &slot) {
        return *this;
      }
      block = slot.block;
      occupied = slot.occupied;
      last_tenant = slot.last_tenant;
      return *this;
    }
    Slot& operator=(Slot&& slot) noexcept
    {
      block = slot.block;
      occupied = slot.occupied;
      last_tenant = slot.last_tenant;
      return *this;
    }

    typename OrderQueue<T>::Block block;
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

  OrderQueue(const size_t queue_size, const size_t block_size)
    : slots(queue_size, Slot(block_size))
    , queue_size(queue_size)
    , block_size(block_size)
  {}

  OrderQueue(const OrderQueue&) = delete;
  OrderQueue(OrderQueue&&) = delete;

protected:
  std::vector<Slot> slots;
  size_t queue_size, block_size;
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
  template<typename T>                                                         \
  class OrderQueue##SUFFIX : public OrderQueue<T>                              \
  {                                                                            \
                                                                               \
  public:                                                                      \
    OrderQueue##SUFFIX(const size_t queue_size, const size_t block_size)       \
      : OrderQueue<T>(queue_size, block_size)                                  \
    {}                                                                         \
                                                                               \
    using Block = typename OrderQueue<T>::Block;                               \
    using Slot = typename OrderQueue<T>::Slot;                                 \
                                                                               \
    void write(Block& block)                                                   \
    {                                                                          \
      PRE_WRITE_LOCK;                                                          \
      const auto num = block.num;                                              \
      auto& target = this->slots[num % this->queue_size];                      \
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
      auto& target = this->slots[this->read_counter % this->queue_size];       \
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
                 &&(num - target.last_tenant <= this->queue_size),
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
                 &&(num - target.last_tenant <= this->queue_size),
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