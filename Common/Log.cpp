#include "Log.h"
#include "Common/Options.h"
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
