#ifndef KMER_ITERATOR_H_
#define KMER_ITERATOR_H_

#include "Common/Sequence.h"
#include "Common/Kmer.h"
#include <limits>
#include <iterator>
#include <string>

struct KmerIterator
: public std::iterator<std::input_iterator_tag, Kmer>
{
	void next()
	{
		for(; m_pos < m_seq.size() - m_k + 1; m_pos++) {
			std::string kmerStr = m_seq.substr(m_pos, m_k);
			size_t pos = kmerStr.find_first_not_of("AGCTagct");
			if (pos != std::string::npos) {
				m_pos += pos;
				continue;
			}
			m_kmer = Kmer(kmerStr);
			if (m_rc)
				m_kmer.reverseComplement();
			return;
		}
		m_pos = std::numeric_limits<std::size_t>::max();
	}

public:

	KmerIterator() :
		m_seq(),
		m_pos(std::numeric_limits<std::size_t>::max()),
		m_kmer() { }

	KmerIterator(const Sequence& seq, unsigned k, bool rc = false)
		: m_seq(seq), m_k(k), m_rc(rc), m_pos(0), m_kmer()
	{
		next();
	}

	const Kmer& operator*() const
	{
		assert(m_pos < m_seq.size() - m_k + 1);
		return m_kmer;
	}

	bool operator==(const KmerIterator& it) const
	{
		return m_pos == it.m_pos;
	}

	bool operator!=(const KmerIterator& it) const
	{
		return !(*this == it);
	}

	KmerIterator& operator++()
	{
		assert(m_pos < m_seq.size() - m_k + 1);
		m_pos++;
		next();
		return *this;
	}

	KmerIterator operator++(int)
	{
		KmerIterator it = *this;
		++*this;
		return it;
	}

	size_t pos() {
		return m_pos;
	}

	static const KmerIterator& end()
	{
		return KmerIterator::m_end;
	}
	
private:

	const Sequence m_seq;
	unsigned m_k;
	bool m_rc;
	size_t m_pos;
	Kmer m_kmer;
	static const KmerIterator m_end;
};

const KmerIterator KmerIterator::m_end = KmerIterator();

#endif
