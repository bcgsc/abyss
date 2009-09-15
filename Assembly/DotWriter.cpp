/** Written by Shaun Jackman <sjackman@bcgsc.ca>. */

#include "DotWriter.h"
#include "AssemblyAlgorithms.h"
#include <ostream>

using namespace std;

/** Return the kmer which are adjacent to this kmer. */
static vector<PackedSeq> getExtensions(
		const PackedSeq& seq, extDirection dir)
{
	vector<PackedSeq> v;
	AssemblyAlgorithms::generateSequencesFromExtension(seq, dir,
			seq.getExtension(dir), v);
	return v;
}

/** Return the kmer which are adjacent to this kmer. */
static vector<PackedSeq> getExtensions(
		const ISequenceCollection& c,
		const PackedSeq& key, extDirection dir)
{
	ExtensionRecord ext;
	int multiplicity = -1;
	c.getSeqData(key, ext, multiplicity);
	vector<PackedSeq> v;
	AssemblyAlgorithms::generateSequencesFromExtension(key, dir,
			ext.dir[dir], v);
	return v;
}

/** Return true if the specified sequence has a single extension in
 * the forward direction, and that extension has a single extension in
 * the reverse direction. If so, return the extension in pSeq.
 */
static bool isContiguous(const ISequenceCollection& c,
		PackedSeq* pSeq, extDirection dir)
{
	vector<PackedSeq> exts = getExtensions(c, *pSeq, dir);
	if (exts.size() != 1)
		return false;
	const PackedSeq& ext = exts[0];
	vector<PackedSeq> rexts = getExtensions(c, ext, !dir);
	if (rexts.size() != 1)
		return false;
	*pSeq = ext;
	return true;
}

static bool isContiguous(const ISequenceCollection& c,
		const PackedSeq& seq, extDirection dir)
{
	PackedSeq ext = seq;
	return isContiguous(c, &ext, dir);
}

/** Find and return the end of the specified contig. Return the length
 * of the contig (counted in kmers) in pLength.
 */
static const PackedSeq findContigEnd(const ISequenceCollection& c,
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
		const ISequenceCollection& c, const PackedSeq& seq)
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
	vector<PackedSeq> exts = getExtensions(seq, SENSE);
	if (exts.size() <= 1)
		return;
	out << seq.decode() << "->{ ";
	for (unsigned i = 0; i < exts.size(); i++)
		out << exts[i].decode() << ' ';
	out << "};\n";
}

/** Write out the contigs that join at the specified sequence. */
static void write_join(ostream& out, const PackedSeq& seq)
{
	vector<PackedSeq> exts = getExtensions(seq, ANTISENSE);
	if (exts.size() <= 1)
		return;
	out << "{ ";
	for (unsigned i = 0; i < exts.size(); i++)
		out << exts[i].decode() << ' ';
	out << "}->" << seq.decode() << ";\n";
}

/** Write out a dot graph around the specified sequence. */
static void write_node(ostream& out,
		const ISequenceCollection& c, const PackedSeq& seq)
{
	write_split(out, seq);
	write_join(out, seq);
	if (!isContiguous(c, seq, ANTISENSE))
		write_contig(out, c, seq);
}

/** Write out a dot graph for the specified collection. */
void DotWriter::write(ostream& out, const ISequenceCollection& c)
{
	out << "digraph g {\n";
	for (ISequenceCollection::const_iterator i = c.begin();
			i != c.end(); ++i) {
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
