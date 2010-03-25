#include "Fcontrol.h"
#include <fcntl.h>

/* Set the FD_CLOEXEC flag of the specified file descriptor. */
int setCloexec(int fd)
{
	int flags = fcntl(fd, F_GETFD, 0);
	if (flags == -1)
		return -1;
	flags |= FD_CLOEXEC;
	return fcntl(fd, F_SETFD, flags);
}
