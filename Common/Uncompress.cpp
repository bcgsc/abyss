/** Uncompress input files using pipes.
 * Hook the standard file opening functions, open, fopen and fopen64.
 * If the extension of the file being opened indicates the file is
 * compressed (.gz, .bz2, .xz), open a pipe to a program that
 * decompresses that file (gunzip, bunzip2 or xzdec) and return a
 * handle to the open pipe.
 * @author Shaun Jackman <sjackman@bcgsc.ca>
 */

#include "config.h"
#if HAVE_LIBDL

#include "Fcntl.h"
#include <cassert>
#include <cstdio> // for perror
#include <cstdlib>
#include <dlfcn.h>
#include <string>
#include <unistd.h>

using namespace std;

/** Tests whether this string ends with the specified suffix. */
static bool endsWith(const string& s, const string& suffix)
{
	ssize_t i = s.length() - suffix.length();
	return i < 0 ? false : s.substr(i) == suffix;
}

static const char* zcatExec(const string& path)
{
	return
		endsWith(path, ".tar") ? "tar -xOf " :
		endsWith(path, ".tar.Z") ? "tar -zxOf " :
		endsWith(path, ".tar.gz") ? "tar -zxOf " :
		endsWith(path, ".tar.bz2") ? "tar -jxOf " :
		endsWith(path, ".tar.xz") ?
			"tar --use-compress-program=xzdec -xOf " :
		endsWith(path, ".Z") ? "gunzip -c" :
		endsWith(path, ".gz") ? "gunzip -c" :
		endsWith(path, ".bz2") ? "bunzip2 -c" :
		endsWith(path, ".xz") ? "xzdec -c" :
		NULL;
}

extern "C" {

/** Open a pipe to uncompress the specified file.
 * @return a file descriptor
 */
static int uncompress(const char *path)
{
	const char *zcat = zcatExec(path);
	assert(zcat != NULL);

	int fd[2];
	if (pipe(fd) == -1)
		return -1;
	int err = setCloexec(fd[0]);
	assert(err == 0);
	(void)err;

	pid_t pid = fork();
	if (pid == -1)
		return -1;

	if (pid == 0) {
		char arg0[16], arg1[16], arg2[16];
		int n = sscanf(zcat, "%s %s %s", arg0, arg1, arg2);
		assert(n == 2 || n == 3);
		dup2(fd[1], STDOUT_FILENO);
		close(fd[1]);
		if (n == 2)
			execlp(arg0, arg0, arg1, path, NULL);
		else
			execlp(arg0, arg0, arg1, arg2, path, NULL);
		perror(zcat);
		exit(EXIT_FAILURE);
	} else {
		close(fd[1]);
		return fd[0];
	}
}

/** Open a pipe to uncompress the specified file.
 * @return a FILE pointer
 */
static FILE* funcompress(const char* path)
{
	int fd = uncompress(path);
	if (fd == -1) {
		perror(path);
		exit(EXIT_FAILURE);
	}
	return fdopen(fd, "r");
}

typedef FILE* (*fopen_t)(const char *path, const char *mode);

/** If the specified file is compressed, return a pipe that
 * uncompresses it.
 */
FILE *fopen(const char *path, const char *mode)
{
	static fopen_t real_fopen;
	if (real_fopen == NULL)
		real_fopen = (fopen_t)dlsym(RTLD_NEXT, "fopen");
	if (real_fopen == NULL) {
		fprintf(stderr, "error: dlsym fopen: %s\n", dlerror());
		exit(EXIT_FAILURE);
	}
	return zcatExec(path) == NULL ? real_fopen(path, mode)
		: funcompress(path);
}

/** If the specified file is compressed, return a pipe that
 * uncompresses it.
 */
FILE *fopen64(const char *path, const char *mode)
{
	static fopen_t real_fopen64;
	if (real_fopen64 == NULL)
		real_fopen64 = (fopen_t)dlsym(RTLD_NEXT, "fopen64");
	if (real_fopen64 == NULL) {
		fprintf(stderr, "error: dlsym fopen64: %s\n", dlerror());
		exit(EXIT_FAILURE);
	}
	return zcatExec(path) == NULL ? real_fopen64(path, mode)
		: funcompress(path);
}

typedef int (*open_t)(const char *path, int flags, mode_t mode);

/** If the specified file is compressed, return a pipe that
 * uncompresses it.
 */
int open(const char *path, int flags, mode_t mode)
{
	static open_t real_open;
	if (real_open == NULL)
		real_open = (open_t)dlsym(RTLD_NEXT, "open");
	if (real_open == NULL) {
		fprintf(stderr, "error: dlsym open: %s\n", dlerror());
		exit(EXIT_FAILURE);
	}
	return zcatExec(path) == NULL ? real_open(path, flags, mode)
		: uncompress(path);
}

} // extern "C"

#endif // HAVE_LIBDL

/** Initialize the uncompress module. */
bool uncompress_init()
{
	return HAVE_LIBDL;
}
