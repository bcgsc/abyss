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
	int pos;
	unsigned mapq;
	std::string cigar;
	std::string mrnm;
	int mpos;
	int isize;
	std::string seq;
	std::string qual;

	/** Flag */
	enum {
		/** the read is paired in sequencing, no matter whether it is
		 * mapped in a pair */
		FPAIRED = 1,
		/** the read is mapped in a proper pair */
		FPROPER_PAIR = 2,
		/** the read itself is unmapped; conflictive with FPROPER_PAIR
		 */
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

	/** Consturct a single-end alignment. */
	SAMRecord(const Alignment& a) :
		qname("*"),
		flag(a.isRC ? FREVERSE : 0),
		rname(a.contig),
		pos(a.contig_start_pos),
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

	/** Construct a paired-end alignment. */
	SAMRecord(const Alignment& a0, const Alignment& a1)
	{
		*this = SAMRecord(a0);
		flag |= FPAIRED;
		if (a1.isRC)
			flag |= FMREVERSE;
		mrnm = a1.contig;
		mpos = a1.contig_start_pos;
		isize = a1.targetAtQueryStart() - a0.targetAtQueryStart();
	}

	/** Parse the specified CIGAR string.
	 * @return an alignment setting the fields read_start_pos,
	 * align_length, and read_length. The other fields will be
	 * uninitialized.
	 */
	static Alignment parseCigar(const std::string& cigar) {
		Alignment a;
		std::istringstream in(cigar);
		unsigned len;
		char type;
		in >> len >> type;
		assert(in.good());
		switch (type) {
			case 'S':
				a.read_start_pos = len;
				in >> len >> type;
				assert(in.good());
				assert(type == 'M');
				a.align_length = len;
				break;
			case 'M':
				a.read_start_pos = 0;
				a.align_length = len;
				break;
			default:
				assert(false);
		}
		if (in >> len >> type) {
			assert(type == 'S');
			(void)in.peek(); // to set the EOF flag
		} else
			len = 0;
		a.read_length = a.read_start_pos + a.align_length + len;
		assert(in.eof());
		return a;
	}

	operator Alignment() const {
		assert(~flag & FUNMAP);
		Alignment a = parseCigar(cigar);
		assert(seq == "*" || (unsigned)a.read_length == seq.length());
		a.contig = rname;
		a.contig_start_pos = pos;
		a.isRC = flag & FREVERSE; // strand of the query
		return a;
	}

	/**
	 * Return the position of the first base of the mate query on the
	 * target extrapolated from the start of the alignment.
	 */
	int mateTargetAtQueryStart() const
	{
		return Alignment(*this).targetAtQueryStart() + isize;
	}

	friend std::ostream& operator <<(std::ostream& out,
			const SAMRecord& o)
	{
		return out << o.qname
			<< '\t' << o.flag
			<< '\t' << o.rname
			<< '\t' << (1 + o.pos)
			<< '\t' << o.mapq
			<< '\t' << o.cigar
			<< '\t' << (o.mrnm == o.rname ? "=" : o.mrnm)
			<< '\t' << (1 + o.mpos)
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
		o.pos--;
		o.mpos--;
		o.qname += o.flag & FREAD1 ? "/1" :
			o.flag & FREAD2 ? "/2" :
			"";
		if (o.mrnm == "=")
			o.mrnm = o.rname;
		return in;
	}

	bool isReverse() const { return flag & FREVERSE; }
	bool isMateReverse() const { return flag & FMREVERSE; }
};

#endif
