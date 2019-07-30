#ifndef STHASH__ITERATOR_H
#define STHASH__ITERATOR_H 1

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

class stHashIterator
{

public:
   
    static std::vector<std::vector<unsigned> > parseSeed(const std::vector<std::string> &seedString) {
        std::vector<std::vector<unsigned> > seedSet;
        for (unsigned i=0; i< seedString.size(); i++) {
            std::vector<unsigned> sSeed;
            for (unsigned j=0; j < seedString[i].size(); j++) {
                if(seedString[i][j]!='1') sSeed.push_back(j);
            }
            seedSet.push_back(sSeed);
        }
        return seedSet;
    }
    
    /**
     * Default constructor. Creates an iterator pointing to
     * the end of the iterator range.
    */
    stHashIterator():
        m_hVec(NULL),
        m_hStn(NULL),
        m_pos(std::numeric_limits<std::size_t>::max())
    {}

    /**
     * Constructor.
     * @param seq address of DNA sequence to be hashed
     * @param seed address of spaced seed
     * @param k k-mer size
     * @param h number of hashes
    */
    stHashIterator(const std::string& seq, const std::vector<std::vector<unsigned> >& seed, unsigned h, unsigned k):
    m_seq(seq), m_seed(seed), m_h(h), m_k(k), m_hVec(new uint64_t[h]), m_hStn(new bool[h]), m_pos(0)
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
        unsigned locN=0;
        while (m_pos<m_seq.length()-m_k+1 && !NTMS64(m_seq.data()+m_pos, m_seed, m_k, m_h, m_fhVal, m_rhVal, locN, m_hVec, m_hStn))
            m_pos+=locN+1;
        if (m_pos >= m_seq.length()-m_k+1)
            m_pos = std::numeric_limits<std::size_t>::max();
    }

    /** Advance iterator right to the next valid k-mer */
    void next()
    {
        ++m_pos;
        if (m_pos >= m_seq.length()-m_k+1) {
            m_pos = std::numeric_limits<std::size_t>::max();
            return;
        }
        if(seedTab[(unsigned char)(m_seq.at(m_pos+m_k-1))]==seedN) {
            m_pos+=m_k;
            init();
        }
        else
            NTMS64(m_seq.data()+m_pos, m_seed, m_seq.at(m_pos-1), m_seq.at(m_pos-1+m_k), m_k, m_h, m_fhVal, m_rhVal, m_hVec, m_hStn);
    }
    
    size_t pos() const{
    	return m_pos;
    }

    /** get pointer to hash values for current k-mer */
    const bool* strandArray() const
    {
        return m_hStn;
    }

    
    /** get pointer to hash values for current k-mer */
    const uint64_t* operator*() const
    {
        return m_hVec;
    }
    

    /** test equality with another iterator */
    bool operator==(const stHashIterator& it) const
    {
        return m_pos == it.m_pos;
    }

    /** test inequality with another iterator */
    bool operator!=(const stHashIterator& it) const
    {
        return !(*this == it);
    }
    
    /** pre-increment operator */
    stHashIterator& operator++()
    {
        next();
        return *this;
    }
    
    /** iterator pointing to one past last element */
    static const stHashIterator end()
    {
        return stHashIterator();
    }
    
    /** destructor */
    ~stHashIterator() {
        if(m_hVec!=NULL) {
            delete [] m_hVec;
            delete [] m_hStn;
        }
    }

private:
    
    /** DNA sequence */
    std::string m_seq;
    
    /** Spaced Seed sequence */
    std::vector<std::vector<unsigned> > m_seed;
   
    /** number of hashes */
    unsigned m_h;

    /** k-mer size */
    unsigned m_k;

    /** hash values */
    uint64_t *m_hVec;

    /** hash strands, forward = 0, reverse-complement = 1 */
    bool *m_hStn;
    
    /** position of current k-mer */
    size_t m_pos;

    /** forward-strand k-mer hash value */
    uint64_t m_fhVal;

    /** reverse-complement k-mer hash value */
    uint64_t m_rhVal;
};

#endif
