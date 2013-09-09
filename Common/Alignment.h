#ifndef ALIGNMENT_H
#define ALIGNMENT_H 1

#include <cassert>
#include <istream>
#include <ostream>
#include <string>

/** An ungapped alignment of a query to a target. */
struct Alignment
{

std::string contig;
int contig_start_pos;
int read_start_pos;
int align_length;
int read_length;
bool isRC;

Alignment() { }

Alignment(const Alignment& o,
		std::string /*qid*/, std::string /*seq*/) { *this = o; }

Alignment(std::string contig, int contig_start, int read_start,
		int align_length, int read_length, bool isRC) :
	contig(contig),
	contig_start_pos(contig_start),
	read_start_pos(read_start),
	align_length(align_length),
	read_length(read_length),
	isRC(isRC)
{
}

/**
 * Return the taret position at the query start.
 * Note: not alignment start, and may be negative
 */
int targetAtQueryStart() const
{
	unsigned tend = contig_start_pos + align_length;
	return !isRC ? contig_start_pos - read_start_pos
		: int(tend + read_start_pos);
}

/** Return the distance between the specified alignments.
 * May be used to calculate fragment size when the alignments are
 * mate pairs.
 */
int operator-(const Alignment& o) const
{
	return targetAtQueryStart() - o.targetAtQueryStart();
}

/** Return an alignment of the reverse complement of the query to
 * the same target.
 */
Alignment flipQuery() const
{
	Alignment rc(*this);
	unsigned qend = read_start_pos + align_length;
	assert(qend <= (unsigned)read_length);
	rc.read_start_pos = read_length - qend;
	rc.isRC = !isRC;
	return rc;
}

static int calculateReverseReadStart(int read_start_pos,
		int read_length, int align_length)
{
	unsigned qend = read_start_pos + align_length;
	return read_length - qend;
}

bool operator<(const Alignment& a1) const
{
	return read_start_pos < a1.read_start_pos;
}

friend std::istream& operator >>(std::istream& in, Alignment& a)
{
	return in >> a.contig
		>> a.contig_start_pos
		>> a.read_start_pos
		>> a.align_length
		>> a.read_length
		>> a.isRC;
}

friend std::ostream& operator <<(std::ostream& out,
		const Alignment& a)
{
	return out << a.contig << ' '
		<< a.contig_start_pos << ' '
		<< a.read_start_pos << ' '
		<< a.align_length << ' '
		<< a.read_length << ' '
		<< a.isRC;
}

}; // struct Alignment

#endif
