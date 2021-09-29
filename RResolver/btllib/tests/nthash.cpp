#include "btllib/nthash.hpp"
#include "helpers.hpp"

#include <iostream>
#include <string>

int
main()
{
  std::string kmer = "ACGTACACTGGACTGAGTCT";

  {
    std::cerr << "Testing single kmer hash values" << std::endl;
    btllib::NtHash nthash(kmer, 3, kmer.size());

    /* Hash values*/
    const std::vector<uint64_t> hashes = { 10434435546371013747U,
                                           16073887395445158014U,
                                           8061578976118370557U };

    nthash.roll();
    TEST_ASSERT_EQ(nthash.get_hash_num(), hashes.size());
    size_t i;
    for (i = 0; i < nthash.get_hash_num(); ++i) {
      TEST_ASSERT_EQ(nthash.hashes()[i], hashes[i]);
    }
    TEST_ASSERT_EQ(i, 3);
  }

  {
    std::cerr << "Testing base substitution" << std::endl;
    btllib::NtHash nthash(kmer, 3, kmer.size());
    std::string kmer_subbed = "ACGCGCACTGGACTGAGTCT";
    btllib::NtHash nthash_subbed(kmer_subbed, 3, kmer_subbed.size());

    nthash.roll();
    nthash.sub({ 3, 4 }, { 'C', 'G' });
    nthash_subbed.roll();
    TEST_ASSERT_EQ(nthash.get_hash_num(), nthash_subbed.get_hash_num());
    size_t i;
    for (i = 0; i < nthash.get_hash_num(); ++i) {
      TEST_ASSERT_EQ(nthash.hashes()[i], nthash_subbed.hashes()[i]);
    }
    TEST_ASSERT_EQ(i, 3);
  }

  {
    std::cerr << "Testing reverse complement" << std::endl;
    /* Reverse complement of kmer*/
    std::string rc_kmer = "AGACTCAGTCCAGTGTACGT";

    btllib::NtHash nthash(kmer, 3, 20);
    btllib::NtHash nthash_rc(rc_kmer, 3, 20);

    nthash.roll();
    nthash_rc.roll();
    TEST_ASSERT_EQ(nthash.get_hash_num(), nthash_rc.get_hash_num());
    size_t i;
    for (i = 0; i < nthash.get_hash_num(); ++i) {
      TEST_ASSERT_EQ(nthash.hashes()[i], nthash_rc.hashes()[i]);
    }
    TEST_ASSERT_EQ(i, 3);
  }

  {
    std::cerr << "Testing rolling hash values" << std::endl;
    btllib::NtHash nthash(kmer, 3, 18);

    /* 18-mers of kmer*/
    std::string kmer1 = "ACGTACACTGGACTGAGT";
    std::string kmer2 = "CGTACACTGGACTGAGTC";
    std::string kmer3 = "GTACACTGGACTGAGTCT";

    std::vector<btllib::NtHash> nthash_vector = {
      btllib::NtHash(kmer1, nthash.get_hash_num(), kmer1.size()),
      btllib::NtHash(kmer2, nthash.get_hash_num(), kmer2.size()),
      btllib::NtHash(kmer3, nthash.get_hash_num(), kmer3.size())
    };

    size_t i;
    for (i = 0; nthash.roll() && nthash_vector[i].roll(); ++i) {
      for (size_t j = 0; j < nthash.get_hash_num(); ++j) {
        TEST_ASSERT_EQ(nthash.hashes()[j], nthash_vector[i].hashes()[j]);
      }
    }
    TEST_ASSERT_EQ(i, 3);
  }

  {
    std::cerr << "Testing spaced seeds" << std::endl;
    std::vector<std::string> seeds = { "11111100000000111111",
                                       "11111111000011111111" };

    /* Point Mutations of Kmer*/
    std::string kmerM1 = "ACGTACACTTGACTGAGTCT";
    std::string kmerM2 = "ACGTACACTGTACTGAGTCT";
    std::string kmerM3 = "ACGTACACTGCACTGAGTCT";
    TEST_ASSERT_EQ(kmerM1.size(), seeds[0].size());
    TEST_ASSERT_EQ(kmerM1.size(), seeds[1].size());

    btllib::SeedNtHash seed_nthash(kmer, seeds, 2, kmer.size());

    std::vector<btllib::SeedNtHash> seed_nthash_vector = {
      btllib::SeedNtHash(kmerM1,
                         seeds,
                         seed_nthash.get_hash_num_per_seed(),
                         seed_nthash.get_k()),
      btllib::SeedNtHash(kmerM2,
                         seeds,
                         seed_nthash.get_hash_num_per_seed(),
                         seed_nthash.get_k()),
      btllib::SeedNtHash(
        kmerM3, seeds, seed_nthash.get_hash_num_per_seed(), seed_nthash.get_k())
    };
    TEST_ASSERT_EQ(seed_nthash.get_hash_num(), seeds.size() * 2);
    TEST_ASSERT_EQ(seed_nthash.get_hash_num(),
                   seed_nthash_vector[0].get_hash_num());

    seed_nthash.roll();
    size_t i;
    for (i = 0; i < seed_nthash_vector.size() && seed_nthash_vector[i].roll();
         i++) {
      for (size_t j = 0; j < seed_nthash.get_hash_num(); j++) {
        TEST_ASSERT_EQ(seed_nthash.hashes()[j],
                       seed_nthash_vector[i].hashes()[j]);
      }
    }
    TEST_ASSERT_EQ(i, 3);
  }

  {
    std::cerr << "Testing RNA" << std::endl;
    btllib::NtHash dna_nthash(kmer, 3, 20);

    std::string rna_kmer = "ACGUACACUGGACUGAGUCU";
    btllib::NtHash rna_nthash(kmer, 3, 20);

    dna_nthash.roll();
    rna_nthash.roll();
    size_t i;
    for (i = 0; i < dna_nthash.get_hash_num(); ++i) {
      TEST_ASSERT_EQ(dna_nthash.hashes()[i], rna_nthash.hashes()[i]);
    }
    TEST_ASSERT_EQ(i, 3);
  }

  return 0;
}