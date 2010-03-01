#ifndef BARRIER_H
#define BARRIER_H 1

#include <unistd.h>
#if _POSIX_BARRIERS > 0

#include <cassert>
#include <cstring> // for memset
#include <pthread.h>

/** Maintain a count of active threads and provide a barrier. */
class Barrier {
  public:
	/** Not thread safe. */
	Barrier(unsigned count)
		: m_threads(count), m_first(false), m_dirty(false)
	{
		pthread_mutex_init(&m_mutex, NULL);
		pthread_rwlock_init(&m_rwlock, NULL);
		assert(count > 0);
		pthread_barrier_init(&m_barrier, NULL, count);
	}

	/** This barrier is not fully constructed. Call init before
	 * calling wait or signal. Not thread safe. */
	Barrier() : m_threads(0), m_first(false), m_dirty(false)
	{
		pthread_mutex_init(&m_mutex, NULL);
		pthread_rwlock_init(&m_rwlock, NULL);
		memset(&m_barrier, 0, sizeof m_barrier);
	}

	/** Destroy this barrier. */
	~Barrier()
	{
		pthread_mutex_destroy(&m_mutex);
		pthread_rwlock_destroy(&m_rwlock);
#if 0
		// The barrier may have already been destroyed by wait.
		pthread_barrier_destroy(&m_barrier);
#endif
	}

	/** Call pthread_barrier_init. Not thread safe. */
	void init(unsigned count)
	{
		m_threads = count;
		m_first = false;
		m_dirty = false;
		assert(count > 0);
		int err = pthread_barrier_init(&m_barrier, NULL, count);
		assert(err == 0);
		(void)err;
	}

	/** Wait for all threads to reach this point. If the number of
	 * threads has changed, it will be reinitialized. */
	void wait()
	{
		pthread_rwlock_rdlock(&m_rwlock);
		pthread_barrier_wait(&m_barrier);
		m_first = true;
		pthread_rwlock_unlock(&m_rwlock);

		pthread_mutex_lock(&m_mutex);
		if (m_first) {
			// Exactly one thread executes this code.
			pthread_rwlock_wrlock(&m_rwlock);
			m_first = false;
			reinit();
			pthread_rwlock_unlock(&m_rwlock);
		}
		pthread_mutex_unlock(&m_mutex);
	}

	/** Create a thread whose sole purpose is to wait on this barrier.
	 */
	void signal(void)
	{
		pthread_t thread;
		pthread_create(&thread, NULL, barrierWait, this);
		pthread_detach(thread);
	}

	/** Set the number of threads. */
	unsigned operator =(unsigned threads)
	{
		pthread_mutex_lock(&m_mutex);
		m_threads = threads;
		m_dirty = true;
		pthread_mutex_unlock(&m_mutex);
		return threads;
	}

	/** Increment the number of threads. */
	unsigned operator ++()
	{
		pthread_mutex_lock(&m_mutex);
		unsigned n = ++m_threads;
		m_dirty = true;
		pthread_mutex_unlock(&m_mutex);
		return n;
	}

	/** Decrement the number of threads. */
	unsigned operator --()
	{
		pthread_mutex_lock(&m_mutex);
		assert(m_threads > 0);
		unsigned n = --m_threads;
		m_dirty = true;
		pthread_mutex_unlock(&m_mutex);
		return n;
	}

  private:
	/** No copy constructor. */
	Barrier(const Barrier&);

	/** Reinitialize this barrier if the number of threads has
	 * changed. Not thread safe.
	 */
	void reinit()
	{
		if (!m_dirty)
			return;
		m_dirty = false;
		int err = pthread_barrier_destroy(&m_barrier);
		assert(err == 0);
		(void)err;
		if (m_threads > 0)
			pthread_barrier_init(&m_barrier, NULL, m_threads);
	}

	/** Wait on the specified barrier. */
	static void* barrierWait(void *arg)
	{
		Barrier* barrier = static_cast<Barrier*>(arg);
		barrier->wait();
		return NULL;
	}

	/** The number of threads associated with this barrier. */
	unsigned m_threads;

	/** This flag is true for the first thread to leave the barrier,
	 * and false for the others. */
	bool m_first;

	/** Whether the barrier needs to be reinitialized. */
	bool m_dirty;

	/** A mutex synchronizing access to this object. */
	pthread_mutex_t m_mutex;

	/** A read-write lock to synchronize wait and destroy. */
	pthread_rwlock_t m_rwlock;

	/** The underlying barrier. */
	pthread_barrier_t m_barrier;
};

#else /* _POSIX_BARRIERS */

class Barrier {
  public:
	void init(unsigned count) { (void)count; }
	void wait() { }
	void signal(void) { }
	unsigned operator =(unsigned threads) { return threads; }
	unsigned operator ++() { return 0; }
	unsigned operator --() { return 0; }
};

#endif /* _POSIX_BARRIERS */

#endif
