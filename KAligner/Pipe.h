#ifndef PIPE_H
#define PIPE_H 1

#include <semaphore.h>
#include <queue>
#include <pthread.h>

/** An asynchronous queue for transmitting data from one thread to
 *  another. */
template <class T>
class Pipe {
  public:
	/** Ready to use after constructed. */
	Pipe(unsigned size = 1024) : open(true)
	{
		assert(size <= SEM_VALUE_MAX);
		assert(size > 0);
		sem_init(&m_sem_in, 0, size);
		sem_init(&m_sem_out, 0, 0);
		pthread_mutex_init(&m_mutex_queue, NULL);
	}

	/** Destoyr semaphores/mutexs. */
	~Pipe()
	{
		sem_destroy(&m_sem_in);
		sem_destroy(&m_sem_out);
		pthread_mutex_destroy(&m_mutex_queue);
	}

	/** Add data to the buffer/queue. */
	void push(T x)
	{
		// Block if pipe is full, or in use.
		sem_wait(&m_sem_in);
		pthread_mutex_lock(&m_mutex_queue);
		assert(open);
		add(x);
		pthread_mutex_unlock(&m_mutex_queue);
		sem_post(&m_sem_out);
	}

	/** Get data and remove it from the buffer. */
	std::pair<T, size_t> pop()
	{
		// block when pipe is empty and open, or in use.
		sem_wait(&m_sem_out);
		pthread_mutex_lock(&m_mutex_queue);

		std::pair<T, size_t> temp = remove();

		pthread_mutex_unlock(&m_mutex_queue);

		// If pipe is open ensure poping will allow one more push.
		// Otherwise, let next accessor know pipe is closed.
		if (temp.second)
			sem_post(&m_sem_in);
		else {
			assert(!open);
			sem_post(&m_sem_out);
		}
		return temp;
	}

	/** Allows a pop when the pipe is empty to signal the pipe is
	 * 	closed. */
	void close()
	{
		pthread_mutex_lock(&m_mutex_queue);
		open = false;
		pthread_mutex_unlock(&m_mutex_queue);
		sem_post(&m_sem_out);
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
	sem_t m_sem_in, m_sem_out;

	/** Mutual exclusion for reading and writing */
	pthread_mutex_t m_mutex_queue;

	/** True if close() has not been called. */
	bool open;

	/** Pipe's buffer */
	std::queue<T> m_queue;
};

#endif
