#ifndef BTLLIB_INDEX_QUEUE_HPP
#define BTLLIB_INDEX_QUEUE_HPP

#include <algorithm>
#include <atomic>
#include <condition_variable>
#include <mutex>
#include <string>
#include <utility>
#include <vector>

namespace btllib {

template<typename T, unsigned QUEUE_SIZE, unsigned BLOCK_SIZE>
class IndexQueue
{

public:
  struct Block
  {

    Block()
      : data(new T[BLOCK_SIZE])
    {}

    Block(const Block&) = delete;

    Block(Block&& block) noexcept
      : current(block.current)
      , count(block.count)
      , index(block.index)
    {
      std::swap(data, block.data);
    }

    Block& operator=(const Block&) = delete;

    Block& operator=(Block&& block) noexcept
    {
      std::swap(data, block.data);
      current = block.current;
      count = block.count;
      index = block.index;
      return *this;
    }

    ~Block() { delete[] data; }

    T* data = nullptr;
    size_t current = 0;
    size_t count = 0;
    size_t index = 0;
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

    typename IndexQueue<T, QUEUE_SIZE, BLOCK_SIZE>::Block block;
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

template<typename T, unsigned QUEUE_SIZE, unsigned BLOCK_SIZE>
class IndexQueueSPMC : public IndexQueue<T, QUEUE_SIZE, BLOCK_SIZE>
{

public:
  using Block = typename IndexQueue<T, QUEUE_SIZE, BLOCK_SIZE>::Block;
  using Slot = typename IndexQueue<T, QUEUE_SIZE, BLOCK_SIZE>::Slot;

  void write(Block& block)
  {
    auto index = block.index;
    auto& target = this->slots[index % QUEUE_SIZE];
    std::unique_lock<std::mutex> busy_lock(target.busy);
    target.occupancy_changed.wait(
      busy_lock, [&] { return !target.occupied || this->closed; });
    if (this->closed) {
      return;
    }
    target.block = std::move(block);
    target.occupied = true;
    target.occupancy_changed.notify_one();
    ++(this->element_count);
  }

  void read(Block& block)
  {
    std::unique_lock<std::mutex> read_lock(read_mutex);

    auto& target = this->slots[this->read_counter % QUEUE_SIZE];
    std::unique_lock<std::mutex> busy_lock(target.busy);
    target.occupancy_changed.wait(
      busy_lock, [&] { return target.occupied || this->closed; });
    if (this->closed) {
      return;
    }
    ++(this->read_counter);

    read_lock.unlock();

    block = std::move(target.block);
    target.occupied = false;
    target.occupancy_changed.notify_one();
    --(this->element_count);
  }

private:
  std::mutex read_mutex;
};

template<typename T, unsigned QUEUE_SIZE, unsigned BLOCK_SIZE>
class IndexQueueSPSC : public IndexQueue<T, QUEUE_SIZE, BLOCK_SIZE>
{

public:
  using Block = typename IndexQueue<T, QUEUE_SIZE, BLOCK_SIZE>::Block;
  using Slot = typename IndexQueue<T, QUEUE_SIZE, BLOCK_SIZE>::Slot;

  void write(Block& block)
  {
    auto index = block.index;
    auto& target = this->slots[index % QUEUE_SIZE];
    std::unique_lock<std::mutex> busy_lock(target.busy);
    target.occupancy_changed.wait(
      busy_lock, [&] { return !target.occupied || this->closed; });
    if (this->closed) {
      return;
    }
    target.block = std::move(block);
    target.occupied = true;
    target.occupancy_changed.notify_one();
    ++(this->element_count);
  }

  void read(Block& block)
  {
    auto& target = this->slots[this->read_counter % QUEUE_SIZE];
    std::unique_lock<std::mutex> busy_lock(target.busy);
    target.occupancy_changed.wait(
      busy_lock, [&] { return target.occupied || this->closed; });
    if (this->closed) {
      return;
    }
    ++(this->read_counter);

    block = std::move(target.block);
    target.occupied = false;
    target.occupancy_changed.notify_one();
    --(this->element_count);
  }
};

template<typename T, unsigned QUEUE_SIZE, unsigned BLOCK_SIZE>
class IndexQueueMPSC : public IndexQueue<T, QUEUE_SIZE, BLOCK_SIZE>
{

public:
  using Block = typename IndexQueue<T, QUEUE_SIZE, BLOCK_SIZE>::Block;
  using Slot = typename IndexQueue<T, QUEUE_SIZE, BLOCK_SIZE>::Slot;

  void write(Block& block)
  {
    auto index = block.index;
    auto& target = this->slots[index % QUEUE_SIZE];
    std::unique_lock<std::mutex> busy_lock(target.busy);
    target.occupancy_changed.wait(busy_lock, [&] {
      return (!target.occupied && (index - target.last_tenant <= QUEUE_SIZE)) ||
             this->closed;
    });
    if (this->closed) {
      return;
    }
    target.block = std::move(block);
    target.occupied = true;
    target.last_tenant = index;
    target.occupancy_changed.notify_all();
    ++(this->element_count);
  }

  void read(Block& block)
  {
    auto& target = this->slots[this->read_counter % QUEUE_SIZE];
    std::unique_lock<std::mutex> busy_lock(target.busy);
    target.occupancy_changed.wait(
      busy_lock, [&] { return target.occupied || this->closed; });
    if (this->closed) {
      return;
    }
    ++(this->read_counter);
    block = std::move(target.block);
    target.occupied = false;
    target.occupancy_changed.notify_all();
    --(this->element_count);
  }
};

} // namespace btllib

#endif