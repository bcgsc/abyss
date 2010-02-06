#ifndef SAM_H
#define SAM_H 1

#include "Aligner.h"
#include <istream>
#include <ostream>
#include <sstream>
#include <string>

struct SAMRecord {
	std::string qname;
	unsigned flag;
	std::string rname;
	unsigned pos;
	unsigned mapq;
	std::string cigar;
	std::string mrnm;
	unsigned mpos;
	int isize;
	std::string seq;
	std::string qual;

	/** Flag */
	enum {
		/** the read is paired in sequencing, no matter whether it is mapped in a pair */
		FPAIRED = 1,
		/** the read is mapped in a proper pair */
		FPROPER_PAIR = 2,
		/** the read itself is unmapped; conflictive with FPROPER_PAIR */
		FUNMAP = 4,
		/** the mate is unmapped */
		FMUNMAP = 8,
		/** the read is mapped to the reverse strand */
		FREVERSE = 16,
		/** the mate is mapped to the reverse strand */
		FMREVERSE = 32,
		/** this is read1 */
		FREAD1 = 64,
		/** this is read2 */
		FREAD2 = 128,
		/** not primary alignment */
		FSECONDARY = 256,
		/** QC failure */
		FQCFAIL = 512,
		/** optical or PCR duplicate */
		FDUP = 1024,
	};

	SAMRecord() { }

	SAMRecord(const Alignment& a) :
		qname("*"),
		flag(a.isRC ? FREVERSE : 0),
		rname(a.contig),
		pos(1 + a.contig_start_pos),
		mapq(255),
		// cigar
		mrnm("*"),
		mpos(0),
		isize(0),
		seq("*"),
		qual("*")
	{
		unsigned qend = a.read_start_pos + a.align_length;
		int qendpad = a.read_length - qend;
		assert(qendpad >= 0);
		std::ostringstream s;
		if (a.read_start_pos > 0)
			s << a.read_start_pos << 'S';
		s << a.align_length << 'M';
		if (qendpad > 0)
			s << qendpad << 'S';
		cigar = s.str();
	}

	operator Alignment() const {
		Alignment a;
		if (flag & FUNMAP)
			return a;

		// Parse the CIGAR string. This could be much improved.
		std::istringstream in(cigar);
		unsigned len;
		char type;
		if (in >> len >> type) {
			switch (type) {
				case 'S':
					a.read_start_pos = len;
					if (in >> len >> type && type == 'M')
						a.align_length = len;
					else
						return a;
					break;
				case 'M':
					a.read_start_pos = 0;
					a.align_length = len;
					break;
				default:
					return a;
			}
		} else
			return a;

		a.contig = rname;
		assert(pos > 0);
		a.contig_start_pos = pos - 1;
		a.read_length = seq.length();
		a.isRC = flag & FREVERSE; // strand of the query
		return a;
	}

	friend std::ostream& operator <<(std::ostream& out,
			const SAMRecord& o)
	{
		return out << o.qname
			<< '\t' << o.flag
			<< '\t' << o.rname
			<< '\t' << o.pos
			<< '\t' << o.mapq
			<< '\t' << o.cigar
			<< '\t' << o.mrnm
			<< '\t' << o.mpos
			<< '\t' << o.isize
			<< '\t' << o.seq
			<< '\t' << o.qual;
	}

	friend std::istream& operator >>(std::istream& in, SAMRecord& o)
	{
		in >> o.qname
			>> o.flag >> o.rname >> o.pos >> o.mapq
			>> o.cigar >> o.mrnm >> o.mpos >> o.isize
			>> o.seq >> o.qual;
		if (!in)
			return in;
		o.qname += o.flag & FREAD1 ? "/1" :
			o.flag & FREAD2 ? "/2" :
			"";
		return in;
	}
};

#endif
