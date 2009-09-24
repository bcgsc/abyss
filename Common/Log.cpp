#include "Log.h"
#include "Common/Options.h"
#include <cstdarg>
#include <iostream>

using namespace std;

/** Print a log message if the verbosity level is at least the
 * specified level.
 */
ostream& logger(int level)
{
	if (opt::verbose < level) {
		static ostream bitBucket(NULL);
		return bitBucket;
	}
	if (opt::rank >= 0)
		cout << opt::rank << ": ";
	return cout;
}

int PrintDebug(int level, const char* format, ...)
{
	if (opt::verbose < level)
		return 0;
	if (opt::rank >= 0)
		printf("%d: ", opt::rank);
	va_list ap;
	va_start(ap, format);
	int retval = vprintf(format, ap);
	va_end(ap);
	return retval;
}
