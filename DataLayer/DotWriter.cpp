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
static bool isContiguous(const ISequenceCollection& c,
		PackedSeq* pSeq)
{
	HitRecord hr = calculateExtension(&c, *pSeq, SENSE);
	if (hr.getNumHits() != 1)
		return false;
	const PackedSeq& ext = hr.getFirstHit().seq;
	HitRecord rhr = calculateExtension(&c, ext, ANTISENSE);
	if (rhr.getNumHits() != 1)
		return false;
	*pSeq = ext;
	return true;
}

/** Find and return the end of the specified contig. Return the length
 * of the contig (counted in kmers) in pLength.
 */
static const PackedSeq findContigEnd(const ISequenceCollection& c,
		const PackedSeq& seq, unsigned* pLength = NULL)
{
	unsigned n = 1;
	PackedSeq cur = seq;
	while (isContiguous(c, &cur))
		n++;
	if (pLength != NULL)
		*pLength = n;
	return cur;
}

/** Write out the specified contig. */
static void write_contig(ostream& out,
		const ISequenceCollection& c, const PackedSeq& seq)
{
	unsigned n;
	out << seq.decode() << "->"
		<< findContigEnd(c, seq, &n).decode()
		<< "[label=\"" << n << "\"];\n";
}

/** Write out the contigs that split at the specified sequence. */
static void write_split(ostream& out,
		const ISequenceCollection& c, const PackedSeq& seq)
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

	for (int i = 0; i < NUM_BASES; i++) {
		char base = BASES[i];
		if (seq.checkExtension(SENSE, base)) {
			ext.setLastBase(SENSE, base);
			write_contig(out, c, ext);
		}
	}
}

/** Write out the contigs that join at the specified sequence. */
static void write_join(ostream& out,
		const ISequenceCollection& c, const PackedSeq& seq)
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
	write_contig(out, c, seq);
}

/** Write out a dot graph for the specified collection. */
void DotWriter::write(ostream& out, const ISequenceCollection& c)
{
	out << "digraph g {\n";
	const SequenceCollectionIterator& end = c.getEndIter();
	for (SequenceCollectionIterator i = c.getStartIter();
			i != end; ++i) {
		if (i->isFlagSet(SF_DELETE))
			continue;
		if (!i->hasExtension(ANTISENSE))
			write_contig(out, c, *i);
		if (i->isAmbiguous(SENSE))
			write_split(out, c, *i);
		if (i->isAmbiguous(ANTISENSE))
			write_join(out, c, *i);
	}
	out << "}" << endl;
}
