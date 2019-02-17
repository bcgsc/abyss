#ifndef MEMORYUTIL_H
#define MEMORYUTIL_H 1

#include "config.h"
#include <cassert>
#include <fstream>
#include <unistd.h> // for sbrk

#if __MACH__
# ifdef __APPLE__
#  include <mach/mach.h> // for mach_task_self
#  include <mach/task.h> // for task_info
# else
extern "C" {
#  include <mach/mach.h> // for mach_task_self and task_info
}
# endif
#endif

/** Return the number of bytes used by the data and stack segments.
 * @return -1 on error
 */
static inline ssize_t getMemoryUsage()
{
#if __MACH__
    struct task_basic_info t_info;
    mach_msg_type_number_t t_info_count = TASK_BASIC_INFO_COUNT;
	int status = task_info(mach_task_self(),
			TASK_BASIC_INFO, (task_info_t)&t_info, &t_info_count);
	assert(status == KERN_SUCCESS);
	return status == KERN_SUCCESS ? (ssize_t)t_info.virtual_size : -1;
#elif HAVE_GETPAGESIZE
	std::ifstream in("/proc/self/statm");
	size_t size, resident, share, text, lib, data;
	return in >> size >> resident >> share >> text >> lib >> data
		? ssize_t(data * getpagesize()) : -1;
#else
	/** Start of the data segment. */
	static intptr_t sbrk0 = reinterpret_cast<intptr_t>(sbrk(0));
	return reinterpret_cast<intptr_t>(sbrk(0)) - sbrk0;
#endif
}

#endif
