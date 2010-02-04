/**
 * Signal handling code, particularly SIGCHLD.
 */

#include "Signal.h"
#include <cassert>
#include <cerrno>
#include <iostream>
#include <pthread.h>
#include <signal.h>
#include <sys/wait.h>

using namespace std;

/** Print the specified exit status. */
static void printStatus(pid_t pid, int status)
{
	if (WIFEXITED(status))
		cerr << "PID " << pid << " exited with status "
			<< WEXITSTATUS(status) << endl;
	else if (WIFSIGNALED(status))
		cerr << "PID " << pid << " killed by signal "
			<< WTERMSIG(status) << endl;
	else
		cerr << "PID " << pid << " exited with code "
			<< status << endl;
}

/** SIGCHLD handler. Reap child processes and report an error if any
 * fail. */
static void sigchldHandler(int sig)
{
	assert(sig == SIGCHLD);
	(void)sig;

	pid_t pid;
	int status;
	while ((pid = waitpid(-1, &status, WNOHANG)) > 0) {
		// Writing to cerr in a signal handler is not allowed, but
		// we're about to exit and an error message would be really
		// helpful.
		if (status != 0) {
			printStatus(pid, status);
			exit(EXIT_FAILURE);
		}
	}
	if (pid == -1 && errno != ECHILD) {
		perror("waitpid");
		exit(EXIT_FAILURE);
	}
}

/** Wait for SIGCHLD signals and call sigchldHandler. */
static void* sigchldThread(void* arg)
{
	sigset_t* sigset = static_cast<sigset_t*>(arg);
	for (;;) {
		int sig;
		int err = sigwait(sigset, &sig);
		assert(err == 0);
		(void)err;
		assert(sig == SIGCHLD);
		sigchldHandler(sig);
	}
	return NULL;
}

/** Start a thread to handle SIGCHLD. */
void signalInit()
{
	static sigset_t sigset;
	sigemptyset(&sigset);
	sigaddset(&sigset, SIGCHLD);
	pthread_sigmask(SIG_BLOCK, &sigset, NULL);

	pthread_t thread;
	pthread_create(&thread, NULL, sigchldThread, &sigset);
	pthread_detach(thread);
}
