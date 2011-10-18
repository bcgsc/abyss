#ifndef SEMAPHORE_H
#define SEMAPHORE_H 1

/** Semaphore class needed since some OS' do not support unnamed
 * semaphores. */
#if __APPLE__
# include <pthread.h>
class Semaphore {
  public:
	Semaphore(unsigned value) : m_value(value)
	{
		pthread_mutex_init(&m_mutex, NULL);
		pthread_cond_init(&m_cv, NULL);
	}

	~Semaphore()
	{
		pthread_mutex_destroy(&m_mutex);
		pthread_cond_destroy(&m_cv);
	}

	void wait()
	{
		pthread_mutex_lock(&m_mutex);
		// Not a spinlock! For some reason signaling a condvar can
		// cause more than one thread waiting to continue...
		while (m_value == 0)
			pthread_cond_wait(&m_cv, &m_mutex);
		assert(m_value > 0);
		m_value--;
		pthread_mutex_unlock(&m_mutex);
	}

	void post()
	{
		pthread_mutex_lock(&m_mutex);
		m_value++;
		pthread_cond_signal(&m_cv);
		pthread_mutex_unlock(&m_mutex);
	}

  private:
	unsigned m_value;
	pthread_mutex_t m_mutex;
	pthread_cond_t m_cv;
};
#else
# include <semaphore.h>
# include <cerrno>
# include <cstring> // for strerror
class Semaphore {
  public:
	Semaphore(unsigned value)
	{
#if SEM_VALUE_MAX
		assert(value <= SEM_VALUE_MAX);
#endif
		if (sem_init(&m_sem, 0, value) == -1) {
			std::cerr << "error: sem_init:" <<
				strerror(errno) << std::endl;
			assert(false);
			exit(EXIT_FAILURE);
		}
	}
	~Semaphore()
	{
		if (sem_destroy(&m_sem) == -1) {
			std::cerr << "error: sem_destroy:" <<
				strerror(errno) << std::endl;
			assert(false);
			exit(EXIT_FAILURE);
		}
	}
	void wait()
	{
		if (sem_wait(&m_sem) == -1) {
			std::cerr << "error: sem_wait:" <<
				strerror(errno) << std::endl;
			assert(false);
			exit(EXIT_FAILURE);
		}
	}
	void post()
	{
		if (sem_post(&m_sem) == -1) {
			std::cerr << "error: sem_post:" <<
				strerror(errno) << std::endl;
			assert(false);
			exit(EXIT_FAILURE);
		}
	}
  private:
	sem_t m_sem;
};
#endif
#endif
