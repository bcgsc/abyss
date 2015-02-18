#ifndef PAIREDDBG_BRANCHRECORD_H
#define PAIREDDBG_BRANCHRECORD_H 1

#include "Common/Options.h" // for opt::verbose

/** Generate the sequence of this contig. */
template <typename It, typename OutIt>
void branchRecordToStr(It it, It last, OutIt out)
{
	assert(it < last);
	std::string k0 = it->first.str();
	std::copy(k0.begin(), k0.end(), out);
	++it;
	Sequence::iterator outa = out + Kmer::length();
	Sequence::iterator outb = out + KmerPair::length();
	for (; it != last; ++it) {
		std::pair<char, char> x = it->first.getLastBaseChar();
		if (*outa == 'N' || *outa == x.first) {
			*outa = x.first;
		} else {
			char amb = ambiguityOr(*outa, x.first);
			if (opt::verbose > 1) {
				std::cerr << "Warning: Expected '" << *outa
					<< "' and saw '" << x.first
					<< "' at " << outa - out << '\n';
				std::copy(out, outb, std::ostream_iterator<char>(std::cerr));
				std::cerr << '\n'
					<< std::string(outa - out - Kmer::length() + 1, ' ')
					<< it->first.str() << '\n'
					<< std::string(outa - out, ' ') << amb << '\n';
			}
			*outa = amb;
		}
		++outa;
		assert(*outb == 'N');
		*outb = x.second;
		++outb;
	}
}

#include "Assembly/BranchRecordBase.h"

#endif
