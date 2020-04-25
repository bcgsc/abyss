/*
 * UnitTests.cpp
 * Unit Tests for various ntHash headers
 *  Created on: January 20, 2020
 *      Author: Johnathan Wong
 */

/* automatically create main() function to run tests */
#define CATCH_CONFIG_MAIN

/* lightweight unit test framework */
#include "vendor/catch.hpp"
#include "ntHashIterator.hpp"
#include "stHashIterator.hpp"

#include <assert.h>
#include <cstring>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <string>

TEST_CASE("test fixture", "[UnitTests]")
{
    /*
	 * NOTES:
	 * - The SECTION blocks below are separate tests that share the
	 * same setup code
	 * - The common setup code is _re-run_ before each SECTION
	 * - In unit test terminology, this type of setup is known as a
	 * "test fixture"
	 * - See
	 * https://github.com/philsquared/Catch/blob/master/docs/tutorial.md#test-cases-and-sections for
	 * details
	 */

    /* START COMMON SETUP CODE */

    std::string kmer = "ACGTACACTGGACTGAGTCT";

    /* END COMMON SETUP CODE */

    SECTION("invariant hash values")
    {
        ntHashIterator invariantIt(kmer, 3, 20);

        /* Hash values*/
        const std::vector<uint64_t> hashes = {10434435546371013747U, 16073887395445158014U, 8061578976118370557U};

        for (int i = 0; i < 3; ++i)
        {
            assert((*invariantIt)[i] == hashes[i]);
        }
    }

    SECTION("reverse complement")
    {
        ntHashIterator kmerIt(kmer, 3, 20);

        /* Reverse complement of kmer*/
        std::string rcKmer = "AGACTCAGTCCAGTGTACGT";
        ntHashIterator rcKmerIt(rcKmer, 3, 20);

        for (int i = 0; i < 3; ++i)
        {
            assert((*kmerIt)[i] == (*rcKmerIt)[i]);
        }
    }

    SECTION("rolling hash values")
    {
        ntHashIterator rollingHashIt(kmer, 3, 18);

        /* 18-mers of kmer*/
        std::string kmer1 = "ACGTACACTGGACTGAGT";
        std::string kmer2 = "CGTACACTGGACTGAGTC";
        std::string kmer3 = "GTACACTGGACTGAGTCT";

        ntHashIterator kmer1It(kmer1, 3, 18);
        ntHashIterator kmer2It(kmer2, 3, 18);
        ntHashIterator kmer3It(kmer3, 3, 18);

        std::vector<ntHashIterator> vecKmerIt = {kmer1It, kmer2It, kmer3It};

        for (int i = 0; i < 3; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                assert((*rollingHashIt)[j] == (*vecKmerIt[i])[j]);
            }
            ++rollingHashIt;
        }
    }

    SECTION("spaced seeds")
    {
        std::vector<std::string> seedString;
        seedString.push_back("11111100000000111111");

        auto numSeeds = seedString.size();
        auto k = seedString[0].size();

        std::vector<std::vector<unsigned>> seedSet = stHashIterator::parseSeed(seedString);

        /* Point Mutations of Kmer*/
        std::string kmerM1 = "ACGTACACTTGACTGAGTCT";
        std::string kmerM2 = "ACGTACACTGTACTGAGTCT";
        std::string kmerM3 = "ACGTACACTGCACTGAGTCT";

        stHashIterator ssItr(kmer, seedSet, numSeeds, 1, k);

        stHashIterator ssM1Itr(kmerM1, seedSet, numSeeds, 1, k);
        stHashIterator ssM2Itr(kmerM2, seedSet, numSeeds, 1, k);
        stHashIterator ssM3Itr(kmerM3, seedSet, numSeeds, 1, k);

        assert((*ssItr)[0] == (*ssM1Itr)[0]);
        assert((*ssItr)[0] == (*ssM2Itr)[0]);
        assert((*ssItr)[0] == (*ssM3Itr)[0]);
    }

    SECTION("RNA")
    {
        ntHashIterator dnaIt(kmer, 3, 20);

        std::string rnaKmer = "ACGUACACUGGACUGAGUCU";
        ntHashIterator rnaIt(kmer, 3, 20);

        for (int i = 0; i < 3; ++i)
        {
            assert((*dnaIt)[i] == (*rnaIt)[i]);
        }
    } /* end test fixture */
}
