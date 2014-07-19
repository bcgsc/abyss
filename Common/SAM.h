#ifndef SAM_H
#define SAM_H 1

#include "config.h" // for SAM_SEQ_QUAL
#include "IOUtil.h"
#include "Alignment.h"
#include "ContigID.h" // for g_contigNames
#include <algorithm> // for swap
#include <cstdlib> // for exit
#include <iostream>
#include <limits> // for numeric_limits
#include <sstream>
#include <string>

namespace opt {
	/** The minimal alignment size. */
	static unsigned minAlign = 1;
}

/** A SAM alignment of a single query. */
struct SAMAlignment {
	std::string rname;
	int pos;
	unsigned short flag;
	unsigned short mapq;
	std::string cigar;

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

	SAMAlignment() :
		rname("*"),
		pos(-1),
		flag(FUNMAP),
		mapq(0) { }

	/** Consturct a single-end alignment. */
	SAMAlignment(const Alignment& a) :
		rname(a.contig),
		pos(a.contig_start_pos),
		flag(a.isRC ? FREVERSE : 0),
		mapq(255)
	{
		unsigned qend = a.read_start_pos + a.align_length;
		int clip0 = a.read_start_pos;
		int clip1 = a.read_length - qend;
		assert(clip1 >= 0);
		if (a.isRC)
			std::swap(clip0, clip1);
		std::ostringstream s;
		if (clip0 > 0)
			s << clip0 << 'S';
		s << a.align_length << 'M';
		if (clip1 > 0)
			s << clip1 << 'S';
		cigar = s.str();
	}

	bool isPaired() const { return flag & FPAIRED; }
	bool isUnmapped() const { return flag & FUNMAP; }
	bool isMateUnmapped() const { return flag & FMUNMAP; }
	bool isReverse() const { return flag & FREVERSE; }
	bool isMateReverse() const { return flag & FMREVERSE; }
	bool isRead1() const { return flag & FREAD1; }
	bool isRead2() const { return flag & FREAD2; }

	/** The alignment coordinates of a gapped alignment. */
	struct CigarCoord {
		/** The length of the query sequence. */
		unsigned qlen;
		/** The start of the alignment on the query. */
		unsigned qstart;
		/** The length of the alignment on the query. */
		unsigned qspan;
		/** The length of the alignment on the target. */
		unsigned tspan;

		/** Parse the specified CIGAR string. */
		CigarCoord(const std::string& cigar)
			: qlen(0), qstart(0), qspan(0), tspan(0)
		{
			if (cigar == "*")
				return;
			std::istringstream in(cigar);
			bool first = true;
			unsigned len;
			char type;
			while (in >> len >> type) {
				switch (type) {
				  case 'H': case 'S':
					if (first)
						qstart = len;
					qlen += len;
					break;
				  case 'M': case 'X': case '=':
					qlen += len;
					qspan += len;
					tspan += len;
					break;
				  case 'I':
					qlen += len;
					qspan += len;
					break;
				  case 'D': case 'N': case 'P':
					tspan += len;
					break;
				  default:
					std::cerr << "error: invalid CIGAR: `"
						<< cigar << "'\n";
					exit(EXIT_FAILURE);
				}
				first = false;
			}
			assert(in.eof());
		}
	};

	/**
	 * Return the position of the first base of the query on the
	 * target extrapolated from the start of the alignment.
	 */
	int targetAtQueryStart() const
	{
		CigarCoord a(cigar);
		assert(a.qstart + a.qspan <= a.qlen);
		return isReverse()
			? pos + a.tspan + (a.qlen - a.qspan - a.qstart)
			: pos - a.qstart;
	}

	/** Parse the specified CIGAR string.
	 * @return an alignment setting the fields read_start_pos,
	 * align_length, and read_length. The other fields will be
	 * uninitialized.
	 */
	static Alignment parseCigar(const std::string& cigar, bool isRC) {
		Alignment a;
		std::istringstream in(cigar);
		unsigned len;
		char type;
		unsigned clip0 = 0;
		a.align_length = 0;
		unsigned qlen = 0;
		unsigned clip1 = 0;
		while (in >> len >> type) {
			switch (type) {
			  case 'I': case 'X': case '=':
				qlen += len;
				clip1 += len;
			  case 'D': case 'N': case 'P':
				if (a.align_length == 0) {
					// Ignore a malformatted CIGAR string whose first
					// non-clipping operation is not M.
					std::cerr << "warning: malformatted CIGAR: "
						<< cigar << std::endl;
				}
				break;
			  case 'M':
				if ((unsigned)a.align_length < len) {
					clip0 += a.align_length + clip1;
					a.align_length = len;
					qlen += len;
					clip1 = 0;
					break;
				}
			  case 'H': case 'S':
				qlen += len;
				clip1 += len;
				break;
			  default:
				std::cerr << "error: invalid CIGAR: `"
					<< cigar << "'\n";
				exit(EXIT_FAILURE);
			}
		}
		a.read_start_pos = isRC ? clip1 : clip0;
		a.read_length = qlen;
		if (!in.eof()){
			std::cerr << "error: invalid CIGAR: `"
				<< cigar << "'\n";
			exit(EXIT_FAILURE);
		}
		return a;
	}

	operator Alignment() const {
		assert(~flag & FUNMAP);
		bool isRC = flag & FREVERSE; // strand of the query
		Alignment a = parseCigar(cigar, isRC);
		a.contig = rname;
		a.contig_start_pos = pos;
		a.isRC = isRC;
		return a;
	}
};

/** A SAM alignment of a query and its mate. */
struct SAMRecord : SAMAlignment {
	std::string qname;
	std::string mrnm;
	int mpos;
	int isize;
#if SAM_SEQ_QUAL
	std::string seq;
	std::string qual;
	std::string tags;
#endif

	/** Consturct a single-end alignment. */
	explicit SAMRecord(const SAMAlignment& a = SAMAlignment(),
			const std::string& qname = "*",
#if SAM_SEQ_QUAL
			const std::string& seq = "*",
			const std::string& qual = "*",
			const std::string& tags = ""
#else
			const std::string& /*seq*/ = "*",
			const std::string& /*qual*/ = "*",
			const std::string& /*tags*/ = ""
#endif
			) :
		SAMAlignment(a),
		qname(qname),
		mrnm("*"),
		mpos(-1),
		isize(0)
#if SAM_SEQ_QUAL
		,
		seq(seq),
		qual(qual),
		tags(tags)
#endif
	{
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

	/** Set the mate mapping fields. */
	void fixMate(const SAMAlignment& o)
	{
		flag &= ~(FPROPER_PAIR | FMUNMAP | FMREVERSE);
		flag |= FPAIRED;
		if (o.isUnmapped())
			flag |= FMUNMAP;
		if (o.isReverse())
			flag |= FMREVERSE;
		mrnm = o.rname;
		mpos = o.pos;
		isize = isMateUnmapped() ? 0
			: o.targetAtQueryStart() - targetAtQueryStart();

		// Fix unaligned mates
		if (!o.isUnmapped() && isUnmapped()) {
			rname = o.rname;
			pos = o.pos;
		} else if (o.isUnmapped() && !isUnmapped()) {
			mrnm = "=";
			mpos = pos;
		}
	}

	void noMate() { flag &= ~FPAIRED; }

	/**
	 * Return the position of the first base of the mate query on the
	 * target extrapolated from the start of the alignment.
	 */
	int mateTargetAtQueryStart() const
	{
		return targetAtQueryStart() + isize;
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
#if SAM_SEQ_QUAL
			<< '\t' << o.seq
			<< '\t' << o.qual;
#else
			<< "\t*\t*";
#endif
	}

	friend std::istream& operator >>(std::istream& in, SAMRecord& o)
	{
		in >> o.qname
			>> o.flag >> o.rname >> o.pos >> o.mapq
			>> o.cigar >> o.mrnm >> o.mpos >> o.isize;
#if SAM_SEQ_QUAL
		in >> o.seq >> o.qual;
		if (in.peek() != '\n')
			in >> o.tags;
#endif
		in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
		if (!in)
			return in;
		o.pos--;
		o.mpos--;
		if (o.mrnm == "=")
			o.mrnm = o.rname;

		// Set the paired flags if qname ends in /1 or /2.
		unsigned l = o.qname.length();
		if (l >= 2 && o.qname[l-2] == '/') {
			switch (o.qname[l-1]) {
				case '1': o.flag |= FPAIRED | FREAD1; break;
				case '2':
				case '3': o.flag |= FPAIRED | FREAD2; break;
				default: return in;
			}
			o.qname.resize(l - 2);
			assert(!o.qname.empty());
		}

		// Set the unmapped flag if the alignment is not long enough.
		CigarCoord a(o.cigar);
		if (a.qspan < opt::minAlign || a.tspan < opt::minAlign)
			o.flag |= FUNMAP;
		return in;
	}
};

/** Set the mate mapping fields of a0 and a1. */
static inline void fixMate(SAMRecord& a0, SAMRecord& a1)
{
	a0.fixMate(a1);
	a1.fixMate(a0);
}

/** Read contig lengths from SAM headers. */
static inline unsigned readContigLengths(std::istream& in, std::vector<unsigned>& lengths)
{
	assert(in);
	assert(lengths.empty());
	assert(g_contigNames.empty());
	for (std::string line; in.peek() == '@' && getline(in, line);) {
		std::istringstream ss(line);
		std::string type;
		ss >> type;
		if (type != "@SQ")
			continue;

		std::string s;
		unsigned len;
		ss >> expect(" SN:") >> s >> expect(" LN:") >> len;
		assert(ss);

		put(g_contigNames, lengths.size(), s);
		lengths.push_back(len);
	}
	if (lengths.empty()) {
		std::cerr << "error: no @SQ records in the SAM header\n";
		exit(EXIT_FAILURE);
	}
	return lengths.size();
}

#endif
