/** Written by Shaun Jackman <sjackman@bcgsc.ca>. */

#include "DotWriter.h"
#include "AssemblyAlgorithms.h"
#include <cassert>
#include <ostream>

using namespace std;

/** Return the kmer which are adjacent to this kmer. */
static vector<Kmer> getExtensions(const SequenceCollectionHash& c,
		const Kmer& key, extDirection dir)
{
	ExtensionRecord ext;
	int multiplicity = -1;
	bool found = c.getSeqData(key, ext, multiplicity);
	assert(found);
	(void)found;
	vector<Kmer> v;
	AssemblyAlgorithms::generateSequencesFromExtension(key, dir,
			ext.dir[dir], v);
	return v;
}

/** Return true if the specified sequence has a single extension in
 * the forward direction, and that extension has a single extension in
 * the reverse direction. If so, return the extension in pKmer.
 */
static bool isContiguous(const SequenceCollectionHash& c,
		Kmer* pKmer, extDirection dir)
{
	vector<Kmer> exts = getExtensions(c, *pKmer, dir);
	if (exts.size() != 1)
		return false;
	const Kmer& ext = exts[0];
	vector<Kmer> rexts = getExtensions(c, ext, !dir);
	if (rexts.size() != 1)
		return false;
	*pKmer = ext;
	return true;
}

static bool isContiguous(const SequenceCollectionHash& c,
		const Kmer& seq, extDirection dir)
{
	Kmer ext = seq;
	return isContiguous(c, &ext, dir);
}

/** Find and return the end of the specified contig. Return the length
 * of the contig (counted in kmers) in pLength.
 */
static const Kmer findContigEnd(const SequenceCollectionHash& c,
		const Kmer& seq, unsigned* pLength = NULL)
{
	unsigned n = 1;
	Kmer cur = seq;
	while (isContiguous(c, &cur, SENSE))
		n++;
	if (pLength != NULL)
		*pLength = n;
	return cur;
}

/** Write out the specified contig. */
static void write_contig(ostream& out,
		const SequenceCollectionHash& c, const Kmer& seq)
{
	unsigned n;
	const Kmer& end = findContigEnd(c, seq, &n);
	if (n == 1)
		return;
	out << seq.str() << "->" << end.str()
		<< "[label=\"" << n << "\"];\n";
}

/** Write out the contigs that split at the specified sequence. */
static void write_split(ostream& out,
		const SequenceCollectionHash& c, const Kmer& seq)
{
	vector<Kmer> exts = getExtensions(c, seq, SENSE);
	if (exts.size() <= 1)
		return;
	out << seq.str() << "->{ ";
	for (unsigned i = 0; i < exts.size(); i++)
		out << exts[i].str() << ' ';
	out << "};\n";
}

/** Write out the contigs that join at the specified sequence. */
static void write_join(ostream& out,
		const SequenceCollectionHash& c, const Kmer& seq)
{
	vector<Kmer> exts = getExtensions(c, seq, ANTISENSE);
	if (exts.size() <= 1)
		return;
	out << "{ ";
	for (unsigned i = 0; i < exts.size(); i++)
		out << exts[i].str() << ' ';
	out << "}->" << seq.str() << ";\n";
}

/** Write out a dot graph around the specified sequence. */
static void write_node(ostream& out,
		const SequenceCollectionHash& c, const Kmer& seq)
{
	write_split(out, c, seq);
	write_join(out, c, seq);
	if (!isContiguous(c, seq, ANTISENSE))
		write_contig(out, c, seq);
}

/** Write out a dot graph for the specified collection. */
void DotWriter::write(ostream& out, const SequenceCollectionHash& c)
{
	out << "digraph g {\n";
	for (SequenceCollectionHash::const_iterator i = c.begin();
			i != c.end(); ++i) {
		if (i->second.deleted())
			continue;
		const Kmer &seq = i->first;
		write_node(out, c, seq);
		Kmer rseq = reverseComplement(seq);
		if (seq != rseq)
			write_node(out, c, rseq);
	}
	out << "}" << endl;
}
