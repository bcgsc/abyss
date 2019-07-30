#ifndef SSHASH__ITERATOR_H
#define SSHASH__ITERATOR_H 1

#include <string>
#include <limits>
#include "nthash.hpp"


/**
 * Iterate over hash values for k-mers in a
 * given DNA sequence.
 *
 * This implementation uses ntHash
 * function to efficiently calculate
 * hash values for successive k-mers.
 */

class ssHashIterator
{

public:

    /**
     * Default constructor. Creates an iterator pointing to
     * the end of the iterator range.
    */
    ssHashIterator():
        m_pos(std::numeric_limits<std::size_t>::max())
    {}

    /**
     * Constructor.
     * @param seq address of DNA sequence to be hashed
     * @param k k-mer size
     * @param h number of hashes
    */
    ssHashIterator(const std::string& seq, const std::vector<bool>& seed, unsigned k):
        m_seq(seq), m_seed(seed), m_k(k), m_hVal(0), m_sVal(0), m_pos(0)
    {
        init();
    }

    /** Initialize internal state of iterator */
    void init()
    {
        if (m_k > m_seq.length()) {
            m_pos = std::numeric_limits<std::size_t>::max();
            return;
        }
        m_sVal = NTS64(m_seq.data()+m_pos, m_seed, m_k, m_hVal);
    }

    /** Advance iterator right to the next valid k-mer */
    void next()
    {
        ++m_pos;
        if (m_pos >= m_seq.length()-m_k+1) {
            m_pos = std::numeric_limits<std::size_t>::max();
            return;
        }
        m_sVal = NTS64(m_seq.data()+m_pos, m_seed, m_seq.at(m_pos-1), m_seq.at(m_pos-1+m_k), m_k, m_hVal);
    }

    size_t pos() const {
        return m_pos;
    }

    /** get pointer to hash values for current k-mer */
    uint64_t operator*() const
    {
        return m_sVal;
    }

    /** test equality with another iterator */
    bool operator==(const ssHashIterator& it) const
    {
        return m_pos == it.m_pos;
    }

    /** test inequality with another iterator */
    bool operator!=(const ssHashIterator& it) const
    {
        return !(*this == it);
    }

    /** pre-increment operator */
    ssHashIterator& operator++()
    {
        next();
        return *this;
    }

    /** iterator pointing to one past last element */
    static const ssHashIterator end()
    {
        return ssHashIterator();
    }

    /** destructor */
    ~ssHashIterator() {
    }

private:

    /** DNA sequence */
    std::string m_seq;

    /** spaced seed */
    std::vector<bool> m_seed;

    /** k-mer size */
    unsigned m_k;

    /** hash value */
    uint64_t m_hVal;

    /** hash value */
    uint64_t m_sVal;

    /** position of current k-mer */
    size_t m_pos;
};

#endif
