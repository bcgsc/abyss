#ifndef GRAPH_OPTIONS_H
#define GRAPH_OPTIONS_H 1

namespace opt {
	/** The size of a k-mer. */
	extern unsigned k;

	/** The file format of the graph when writing. */
	extern int format;
}

/** Enumeration of output formats */
enum { ADJ, ASQG, DIST, DOT, DOT_MEANCOV, GFA, SAM, TSV };

#endif
