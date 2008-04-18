/** Written by Shaun Jackman <sjackman@bcgsc.ca>. */

#include "DotWriter.h"
#include "AssemblyAlgorithms.h"
#include "CommonDefs.h"
#include <ostream>

using namespace std;

/** Return true if the specified sequence has a single extension in
 * the forward direction, and that extension has a single extension in
 * the reverse direction. If so, return the extension in pSeq.
 */
static bool isContiguous(ISequenceCollection& c,
		PackedSeq* pSeq, extDirection dir)
{
	HitRecord hr = calculateExtension(&c, *pSeq, dir);
	if (hr.getNumHits() != 1)
		return false;
	const PackedSeq& ext = hr.getFirstHit().seq;
	HitRecord rhr = calculateExtension(&c, ext, !dir);
	if (rhr.getNumHits() != 1)
		return false;
	*pSeq = ext;
	return true;
}

static bool isContiguous(ISequenceCollection& c,
		const PackedSeq& seq, extDirection dir)
{
	PackedSeq ext = seq;
	return isContiguous(c, &ext, dir);
}

/** Find and return the end of the specified contig. Return the length
 * of the contig (counted in kmers) in pLength.
 */
static const PackedSeq findContigEnd(ISequenceCollection& c,
		const PackedSeq& seq, unsigned* pLength = NULL)
{
	unsigned n = 1;
	PackedSeq cur = seq;
	while (isContiguous(c, &cur, SENSE))
		n++;
	if (pLength != NULL)
		*pLength = n;
	return cur;
}

/** Write out the specified contig. */
static void write_contig(ostream& out,
		ISequenceCollection& c, const PackedSeq& seq)
{
	unsigned n;
	const PackedSeq& end = findContigEnd(c, seq, &n);
	if (n == 1)
		return;
	out << seq.decode() << "->" << end.decode()
		<< "[label=\"" << n << "\"];\n";
}

/** Write out the contigs that split at the specified sequence. */
static void write_split(ostream& out, const PackedSeq& seq)
{
	out << seq.decode() << "->{ ";
	PackedSeq ext(seq);
	ext.rotate(SENSE, 'A');
	for (int i = 0; i < NUM_BASES; i++) {
		char base = BASES[i];
		if (seq.checkExtension(SENSE, base)) {
			ext.setLastBase(SENSE, base);
			out << ext.decode() << " ";
		}
	}
	out << "};\n";
}

/** Write out the contigs that join at the specified sequence. */
static void write_join(ostream& out, const PackedSeq& seq)
{
	out << "{ ";
	PackedSeq ext(seq);
	ext.rotate(ANTISENSE, 'A');
	for (int i = 0; i < NUM_BASES; i++) {
		char base = BASES[i];
		if (seq.checkExtension(ANTISENSE, base)) {
			ext.setLastBase(ANTISENSE, base);
			out << ext.decode() << " ";
		}
	}
	out << "}->" << seq.decode() << ";\n";
}

/** Write out a dot graph for the specified collection. */
void DotWriter::write(ostream& out, ISequenceCollection& c)
{
	out << "digraph g {\n";
	const SequenceCollectionIterator& end = c.getEndIter();
	for (SequenceCollectionIterator i = c.getStartIter();
			i != end; ++i) {
		if (i->isFlagSet(SF_DELETE))
			continue;
		if (i->isAmbiguous(SENSE))
			write_split(out, *i);
		if (i->isAmbiguous(ANTISENSE))
			write_join(out, *i);
		if (!isContiguous(c, *i, ANTISENSE))
			write_contig(out, c, *i);
	}
	out << "}" << endl;
}
