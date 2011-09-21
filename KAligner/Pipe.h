#ifndef PIPE_H
#define PIPE_H 1

#include "Semaphore.h"
#include <queue>
#include <pthread.h>

/** An asynchronous queue for transmitting data from one thread to
 *  another. */
template <class T>
class Pipe {
  public:
	/** Ready to use after constructed. Not thread safe. */
	Pipe(unsigned size = 1024)
		: m_sem_in(size), m_sem_out(0), m_open(true)
	{
		assert(size > 0);
		pthread_mutex_init(&m_mutex_queue, NULL);
	}

	/** Destoyr semaphores/mutexs. Not thread safe. */
	~Pipe()
	{
		pthread_mutex_destroy(&m_mutex_queue);
	}

	/** Add data to the buffer/queue. */
	void push(T x)
	{
		// Block if pipe is full, or in use.
		m_sem_in.wait();
		pthread_mutex_lock(&m_mutex_queue);
		assert(m_open);
		add(x);
		pthread_mutex_unlock(&m_mutex_queue);
		m_sem_out.post();
	}

	/** Get data and remove it from the buffer. */
	std::pair<T, size_t> pop()
	{
		// block when pipe is empty and m_open, or in use.
		m_sem_out.wait();
		pthread_mutex_lock(&m_mutex_queue);

		std::pair<T, size_t> temp = remove();

		pthread_mutex_unlock(&m_mutex_queue);

		// If pipe is m_open ensure poping will allow one more push.
		// Otherwise, let next accessor know pipe is closed.
		if (temp.second)
			m_sem_in.post();
		else {
			assert(!m_open);
			m_sem_out.post();
		}
		return temp;
	}

	/** Allows a pop when the pipe is empty to signal the pipe is
	 * 	closed. */
	void close()
	{
		pthread_mutex_lock(&m_mutex_queue);
		m_open = false;
		pthread_mutex_unlock(&m_mutex_queue);
		m_sem_out.post();
	}

  private:
  	/** Add an element to the buffer. */
  	void add(const T& t) { m_queue.push(t); }

  	/** Remove an element from the buffer. */
	std::pair<T, size_t> remove()
	{
		std::pair<T, size_t> temp;
		if (!m_queue.empty()) {
			temp.first = m_queue.front();
			temp.second = 1;
			m_queue.pop();
		} else {
			temp.second = 0;
		}
		return temp;
	}

	/** Semaphores to block read on empty, or write on full. */
	Semaphore m_sem_in, m_sem_out;

	/** Mutual exclusion for reading and writing */
	pthread_mutex_t m_mutex_queue;

	/** True if close() has not been called. */
	bool m_open;

	/** Pipe's buffer */
	std::queue<T> m_queue;
};

#endif
