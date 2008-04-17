#include "DotWriter.h"
#include "CommonDefs.h"
#include <algorithm>
#include <ostream>

using namespace std;

static const PackedSeq& findSeq(const ISequenceCollection &c,
		const PackedSeq &seq)
{
	SequenceCollectionIterator i = find(
			c.getStartIter(), c.getEndIter(),
			seq);
	return *i;
}

static const PackedSeq& getFirstExt(const ISequenceCollection &c,
		const PackedSeq &seq, extDirection dir)
{
	for (int i = 0; i < NUM_BASES; i++) {
		char base = BASES[i];
		if (seq.checkExtension(dir, base)) {
			PackedSeq ext(seq);
			ext.rotate(dir, base);
			return findSeq(c, ext);
		}
	}
	assert(false);
}

static const PackedSeq& findContigEnd(const ISequenceCollection &c,
		const PackedSeq &seq, extDirection dir,
		unsigned *pn = NULL)
{
	unsigned n = 1;
	const PackedSeq *p = &seq;
	while (p->hasExtension(dir) && !p->isAmbiguous(dir)) {
		const PackedSeq *ext = &getFirstExt(c, *p, dir);
		extDirection rdir = dir == SENSE ? ANTISENSE : SENSE;
		if (ext->isAmbiguous(rdir))
			break;
		n++;
		p = ext;
	}
	if (pn != NULL)
		*pn = n;
	return *p;
}

static void write_contig(ostream &out,
		const ISequenceCollection &c, const PackedSeq &seq)
{
	unsigned n;
	out << seq.decode() << "->"
		<< findContigEnd(c, seq, SENSE, &n).decode()
		<< "[label=\"" << n << "\"];\n";
}

static void write_split(ostream &out,
		const ISequenceCollection &c, const PackedSeq &seq)
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
			write_contig(out, c, findSeq(c, ext));
		}
	}
}

static void write_join(ostream &out,
		const ISequenceCollection &c, const PackedSeq &seq)
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

void DotWriter::write(ostream &out, const ISequenceCollection& c)
{
	out << "digraph g {\n";
	for (SequenceCollectionIterator i = c.getStartIter(),
			end = c.getEndIter(); i != end; ++i) {
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
