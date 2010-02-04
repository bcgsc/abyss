/**
 * Signal handling code, particularly SIGCHLD.
 */

#include "Signal.h"
#include <cassert>
#include <iostream>
#include <signal.h>
#include <sys/wait.h>

using namespace std;

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

/** Install a handler for SIGCHLD. */
void signalInit()
{
	struct sigaction action;
	action.sa_handler = sigchldHandler;
	sigemptyset(&action.sa_mask);
	action.sa_flags = SA_RESTART;
	sigaction(SIGCHLD, &action, NULL);
}
