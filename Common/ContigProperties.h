#ifndef CONTIGPROPERTIES_H
#define CONTIGPROPERTIES_H 1

#include <cassert>
#include <istream>
#include <ostream>

namespace opt {
	extern int k;

	/** Output format. */
	extern int format;
}

/** Enumeration of output formats */
enum { ADJ, DOT, SAM };

/** Contig properties. */
struct ContigProperties {
	unsigned length;
	unsigned coverage;

	ContigProperties() : length(0), coverage(0) { }
	ContigProperties(unsigned length, unsigned coverage)
		: length(length), coverage(coverage) { }

	ContigProperties& operator +=(const ContigProperties& o)
	{
		length += o.length - opt::k + 1;
		coverage += o.coverage;
		return *this;
	}

	ContigProperties operator +(const ContigProperties& o) const
	{
		return ContigProperties(*this) += o;
	}

	friend std::ostream& operator <<(std::ostream& out,
			const ContigProperties& o)
	{
		float coverage = opt::k <= 0 ? 0
			: (float)o.coverage / (o.length - opt::k + 1);
		switch (opt::format) {
		  case ADJ:
			return out << ' ' << o.length << ' ' << o.coverage;
		  case DOT:
			out << "l=" << o.length;
			return opt::k > 0
				? (out << " c=" << coverage)
				: (out << " C=" << o.coverage);
		  case SAM:
			out << "\tLN:" << o.length;
			return opt::k > 0
				? (out << "\tXc:" << coverage)
				: (out << "\tXC:" << o.coverage);
		}
		return out;
	}

	friend std::istream& operator >>(std::istream& in,
			ContigProperties& o)
	{
		return in >> o.length >> o.coverage;
	}
};

/** The distance between two contigs. */
struct Distance {
	int distance;
	Distance() : distance(-opt::k + 1) { }
	Distance(int d) : distance(d) { }

	bool operator==(const Distance& o) const
	{
		return distance == o.distance;
	}

	bool operator!=(const Distance& o) const
	{
		return distance != o.distance;
	}

	friend std::ostream& operator<<(std::ostream& out,
			const Distance& o)
	{
		return out << "d=" << o.distance;
	}

	friend std::istream& operator>>(std::istream& in, Distance& o)
	{
		char c0, c1;
		in >> c0 >> c1;
		assert(c0 == 'd');
		assert(c1 == '=');
		return in >> o.distance;
	}
};

#endif
