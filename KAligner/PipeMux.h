#ifndef PIPEMUX_H
#define PIPEMUX_H 1

#include "Pipe.h"
#include <semaphore.h>
#include <pthread.h>
#include <vector>

template <class T>
class PipeMux {
  public:
  	/** Default constructor. */
  	PipeMux(size_t pipe_size = 1)
		: m_index(0), m_entry_num(0), m_pipe_size(pipe_size)
	{
		pthread_rwlock_init(&m_rwlock_vecs, NULL);
		pthread_mutex_init(&m_mutex_index, NULL);
	}

	/** Destroy remaining pipes, and mutexes. */
	~PipeMux()
	{
		pthread_rwlock_destroy(&m_rwlock_vecs);
		pthread_mutex_destroy(&m_mutex_index);
		assert(m_pipes.empty());
		assert(m_mutex_pipes.empty());
	}

	/** Instantiates a new pipe and adds it to this PipeMux. */
	Pipe<T>* addPipe()
	{
		Pipe<T>* p = new Pipe<T>(m_pipe_size);
		pthread_mutex_t* m = new pthread_mutex_t;
		pthread_mutex_init(m, NULL);

		pthread_rwlock_wrlock(&m_rwlock_vecs);
		
		m_pipes.push_back(p);
		m_mutex_pipes.push_back(m);
		
		pthread_rwlock_unlock(&m_rwlock_vecs);
		
		return p;
	}

	/** Returns the next value from the appropriate pipe, deletes
	 * closed pipes and returns */
	std::pair<T, size_t> nextValue()
	{
		size_t entry;
		std::pair<T, size_t> t(T(), 0);
		do {
			pthread_rwlock_rdlock(&m_rwlock_vecs);
			pthread_mutex_lock(&m_mutex_index);

			if (m_pipes.empty() && m_mutex_pipes.empty()) {
				pthread_rwlock_unlock(&m_rwlock_vecs);
				pthread_mutex_unlock(&m_mutex_index);
				return t;
			}
			unsigned i = m_index;
			m_index = m_index + 1 < m_pipes.size() ? m_index + 1 : 0;
			entry = ++m_entry_num;

			assert(i < m_mutex_pipes.size());
			pthread_mutex_lock(m_mutex_pipes[i]);
			pthread_mutex_unlock(&m_mutex_index);

			assert(i < m_pipes.size());
			Pipe<T>* p_pipe = m_pipes[i];
			t = p_pipe->pop();

			// you know you're fed up with race conditions when...
			assert(i < m_mutex_pipes.size());
			pthread_mutex_unlock(m_mutex_pipes[i]);
			pthread_rwlock_unlock(&m_rwlock_vecs);
			if (!t.second)
				removePipe(p_pipe, entry);
		} while (!t.second);
		t.second = entry;
		return t;
	}

	bool invalidEntry(size_t e)
	{
		pthread_rwlock_rdlock(&m_rwlock_vecs);
		for (unsigned i = 0; i < m_invalid_entries.size(); i++) {
			assert(i < m_invalid_entries.size());
			if (m_invalid_entries[i] == e) {
				pthread_rwlock_unlock(&m_rwlock_vecs);
				return true;
			}
		}
		pthread_rwlock_unlock(&m_rwlock_vecs);
		return false;
	}

	/** Checks that the PipeMux is empty. */
	bool empty() {
		pthread_rwlock_rdlock(&m_rwlock_vecs);
		bool isEmpty = m_pipes.empty();
		pthread_rwlock_unlock(&m_rwlock_vecs);
		return isEmpty;
	}

  private:
	std::vector<Pipe<T>*> m_pipes;
	std::vector<pthread_mutex_t*> m_mutex_pipes;
	std::vector<size_t> m_invalid_entries;
	pthread_rwlock_t m_rwlock_vecs;
	pthread_mutex_t m_mutex_index;
	unsigned m_index;
	size_t m_entry_num;
	size_t m_pipe_size;

	/** Removes Pipe p if it is still present in m_pipes. */
	void removePipe(Pipe<T>* p, size_t entry)
	{
		pthread_rwlock_wrlock(&m_rwlock_vecs);
		m_invalid_entries.push_back(entry);
		unsigned i;
		for (i = 0; i < m_pipes.size(); i++) {
			assert(i < m_pipes.size());
			if (m_pipes[i] == p)
				break;
		}
		if (i >= m_pipes.size()) {
			pthread_rwlock_unlock(&m_rwlock_vecs);
			return;
		}
		assert(i < m_pipes.size());
		delete m_pipes[i];
		m_pipes.erase(m_pipes.begin()+i);
		assert(i < m_mutex_pipes.size());
		pthread_mutex_destroy(m_mutex_pipes[i]);
		delete m_mutex_pipes[i];
		m_mutex_pipes.erase(m_mutex_pipes.begin()+i);
		// Make sure the index points to the next element.
		pthread_mutex_lock(&m_mutex_index);
		m_index = m_index == m_pipes.size() ? 0 : m_index;
		pthread_mutex_unlock(&m_mutex_index);

		pthread_rwlock_unlock(&m_rwlock_vecs);
	}

};

#endif
