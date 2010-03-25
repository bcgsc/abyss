/**
 * Signal handling code, particularly SIGCHLD.
 */

#include "SignalHandler.h"
#include <cassert>
#include <cerrno>
#include <cstdio> // for perror
#include <cstdlib>
#include <iostream>
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

/** Install a handler for SIGCHLD. */
void signalInit()
{
	struct sigaction action;
	action.sa_handler = sigchldHandler;
	sigemptyset(&action.sa_mask);
	action.sa_flags = SA_RESTART;
	sigaction(SIGCHLD, &action, NULL);
}
