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

#include "Fcontrol.h"
#include "SignalHandler.h"
#include "StringUtil.h"
#include <cassert>
#include <cstdio> // for perror
#include <cstdlib>
#include <dlfcn.h>
#include <string>
#include <unistd.h>

using namespace std;

static const char* wgetExec(const string& path)
{
	return
		startsWith(path, "http://") ? "wget -O-" :
		startsWith(path, "https://") ? "wget -O-" :
		startsWith(path, "ftp://") ? "wget -O-" :
		NULL;
}

static const char* zcatExec(const string& path)
{
	return
		endsWith(path, ".ar") ? "ar -p" :
		endsWith(path, ".tar") ? "tar -xOf" :
		endsWith(path, ".tar.Z") ? "tar -zxOf" :
		endsWith(path, ".tar.gz") ? "tar -zxOf" :
		endsWith(path, ".tar.bz2") ? "tar -jxOf" :
		endsWith(path, ".tar.xz") ?
			"tar --use-compress-program=xzdec -xOf" :
		endsWith(path, ".Z") ? "gunzip -c" :
		endsWith(path, ".gz") ? "gunzip -c" :
		endsWith(path, ".bz2") ? "bunzip2 -c" :
		endsWith(path, ".xz") ? "xzdec -c" :
		endsWith(path, ".zip") ? "unzip -p" :
		endsWith(path, ".bam") ? "samtools view -h" :
		endsWith(path, ".jf") ? "jellyfish dump" :
		endsWith(path, ".jfq") ? "jellyfish qdump" :
		endsWith(path, ".sra") ? "fastq-dump -Z --split-spot" :
		endsWith(path, ".url") ? "wget -O- -i" :
		NULL;
}

extern "C" {

/** Open a pipe to uncompress the specified file.
 * Not thread safe.
 * @return a file descriptor
 */
static int uncompress(const char *path)
{
	const char *wget = wgetExec(path);
	const char *zcat = wget != NULL ? wget : zcatExec(path);
	assert(zcat != NULL);

	int fd[2];
	if (pipe(fd) == -1)
		return -1;
	int err = setCloexec(fd[0]);
	assert(err == 0);
	(void)err;

	char arg0[16], arg1[16], arg2[16];
	int n = sscanf(zcat, "%s %s %s", arg0, arg1, arg2);
	assert(n == 2 || n == 3);

	/* It would be more portable to use fork than vfork, but fork can
	 * fail with ENOMEM when the process calling fork is using a lot
	 * of memory. A workaround for this problem is to set
	 * sysctl vm.overcommit_memory=1
	 */
#if HAVE_WORKING_VFORK
	pid_t pid = vfork();
#else
	pid_t pid = fork();
#endif
	if (pid == -1)
		return -1;

	if (pid == 0) {
		dup2(fd[1], STDOUT_FILENO);
		close(fd[1]);
		if (n == 2)
			execlp(arg0, arg0, arg1, path, NULL);
		else
			execlp(arg0, arg0, arg1, arg2, path, NULL);
		// Calling perror after vfork is not allowed, but we're about
		// to exit and an error message would be really helpful.
		perror(arg0);
		_exit(EXIT_FAILURE);
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

	// open a web address
	if (wgetExec(path) != NULL)
		return funcompress(path);
	
	// to check if the file exists, we need to attempt to open it
	FILE* stream = real_fopen(path, mode);
	if (!stream || zcatExec(path) == NULL)
		return stream;
	else {
		fclose(stream);
		return funcompress(path);
	}
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

	// open a web address
	if (wgetExec(path) != NULL)
		return funcompress(path);
	
	// to check if the file exists, we need to attempt to open it
	FILE* stream = real_fopen64(path, mode);
	if (!stream || zcatExec(path) == NULL)
		return stream;
	else {
		fclose(stream);
		return funcompress(path);
	}
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

	// open a web address
	if (wgetExec(path) != NULL)
		return uncompress(path);
	
	// to check if the file exists, we need to attempt to open it
	int filedesc = real_open(path, flags, mode);
	if (filedesc < 0 || zcatExec(path) == NULL)
		return filedesc;
	else {
		close(filedesc);
		return uncompress(path);
	}
}

} // extern "C"

#endif // HAVE_LIBDL

/** Initialize the uncompress module. */
bool uncompress_init()
{
#if HAVE_LIBDL
	signalInit();
#endif
	return HAVE_LIBDL;
}
