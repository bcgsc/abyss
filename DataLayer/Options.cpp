#include "DataLayer/Options.h"

namespace opt {
	/** Discard reads that failed the chastity filter. */
	int chastityFilter = 1;

	/** Trim masked (lower case) characters from the ends of
	 * sequences.
	 */
	int trimMasked = 1;
}
