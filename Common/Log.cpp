#include "Log.h"
#include "Options.h"
#include <cstdarg>

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
