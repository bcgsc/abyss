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

#include <cassert>
#include <cstdio> // for perror
#include <cstdlib>
#include <dlfcn.h>
#include <iostream>
#include <signal.h>
#include <string>
#include <sys/wait.h>
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

static string zcatCommand(const string& path)
{
	const char *zcat = zcatExec(path);
	return zcat == NULL ? "" : string(zcat) + ' ' + path;
}

extern "C" {

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

	string command = zcatCommand(string(path));
	return command.empty() ? real_fopen(path, mode)
		: popen(command.c_str(), mode);
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

	string command = zcatCommand(string(path));
	return command.empty() ? real_fopen64(path, mode)
		: popen(command.c_str(), mode);
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

	const char *zcat = zcatExec(path);
	if (zcat == NULL)
		return real_open(path, flags, mode);

	int fd[2];
	if (pipe(fd) == -1)
		return -1;

	pid_t pid = fork();
	if (pid == -1)
		return -1;

	if (pid == 0) {
		char arg0[16], arg1[16], arg2[16];
		int n = sscanf(zcat, "%s %s %s", arg0, arg1, arg2);
		assert(n == 2 || n == 3);
		close(fd[0]);
		dup2(fd[1], STDOUT_FILENO);
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

} // extern "C"

/** SIGCHLD handler. Reap child processes and report an error if any
 * fail. */
static void sigchldHandler(int sig)
{
	assert(sig == SIGCHLD);
	(void)sig;
	int status;
	pid_t pid = wait(&status);
	if (pid == -1) {
		perror("waitpid");
		exit(EXIT_FAILURE);
	}
	// Writing to cerr in a signal handler is not allowed, but we're
	// about to exit and an error message would be really helpful.
	if (status != 0) {
		if (WIFEXITED(status))
			cerr << "PID " << pid << " exited with status "
				<< WEXITSTATUS(status) << endl;
		else if (WIFSIGNALED(status))
			cerr << "PID " << pid << " killed by signal "
				<< WTERMSIG(status) << endl;
		else
			cerr << "PID " << pid << " exited with code "
				<< status << endl;
		exit(EXIT_FAILURE);
	}
}

#endif // HAVE_LIBDL

/** Initialize the uncompress module. Install a SIGCHLD handler. */
bool uncompress_init()
{
#if HAVE_LIBDL
	cerr << "SIGCHLD handler installed" << endl;
	signal(SIGCHLD, sigchldHandler);
	return true;
#else
	return false;
#endif
}
