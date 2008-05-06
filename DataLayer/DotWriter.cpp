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
	HitRecord hr = AssemblyAlgorithms::calculateExtension(&c, *pSeq, dir);
	if (hr.getNumHits() != 1)
		return false;
	const PackedSeq& ext = hr.getFirstHit().seq;
	HitRecord rhr = AssemblyAlgorithms::calculateExtension(&c, ext, !dir);
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
static void write_split(ostream& out,
		ISequenceCollection& c, const PackedSeq& seq)
{
	HitRecord hr = AssemblyAlgorithms::calculateExtension(&c, seq, SENSE);
	unsigned hits = hr.getNumHits();
	if (hits <= 1)
		return;
	out << seq.decode() << "->{ ";
	for (unsigned i = 0; i < hits; i++)
		out << hr.getHit(i).seq.decode() << ' ';
	out << "};\n";
}

/** Write out the contigs that join at the specified sequence. */
static void write_join(ostream& out,
		ISequenceCollection& c, const PackedSeq& seq)
{
	HitRecord hr = AssemblyAlgorithms::calculateExtension(&c, seq, ANTISENSE);
	unsigned hits = hr.getNumHits();
	if (hits <= 1)
		return;
	out << "{ ";
	for (unsigned i = 0; i < hits; i++)
		out << hr.getHit(i).seq.decode() << ' ';
	out << "}->" << seq.decode() << ";\n";
}

/** Write out a dot graph around the specified sequence. */
static void write_node(ostream& out,
		ISequenceCollection& c, const PackedSeq& seq)
{
	write_split(out, c, seq);
	write_join(out, c, seq);
	if (!isContiguous(c, seq, ANTISENSE))
		write_contig(out, c, seq);
}

/** Write out a dot graph for the specified collection. */
void DotWriter::write(ostream& out, ISequenceCollection& c)
{
	out << "digraph g {\n";
	const SequenceCollectionIterator& end = c.getEndIter();
	for (SequenceCollectionIterator i = c.getStartIter();
			i != end; ++i) {
		const PackedSeq &seq = *i;
		if (seq.isFlagSet(SF_DELETE))
			continue;
		write_node(out, c, seq);
		PackedSeq rseq = reverseComplement(seq);
		if (seq != rseq)
			write_node(out, c, rseq);
	}
	out << "}" << endl;
}
