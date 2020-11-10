#include "btllib/nthash.hpp"

#include <cassert>
#include <iostream>
#include <string>

int
main()
{
  std::string kmer = "ACGTACACTGGACTGAGTCT";

  {
    std::cerr << "Testing single kmer hash values" << std::endl;
    btllib::NtHash nthash(kmer, kmer.size(), 3);

    /* Hash values*/
    const std::vector<uint64_t> hashes = { 10434435546371013747U,
                                           16073887395445158014U,
                                           8061578976118370557U };

    nthash.roll();
    assert(nthash.get_hash_num() == hashes.size());
    size_t i;
    for (i = 0; i < nthash.get_hash_num(); ++i) {
      assert(nthash.hashes()[i] == hashes[i]);
    }
    assert(i == 3);
  }

  {
    std::cerr << "Testing base substitution" << std::endl;
    btllib::NtHash nthash(kmer, kmer.size(), 3);
    std::string kmer_subbed = "ACGCGCACTGGACTGAGTCT";
    btllib::NtHash nthash_subbed(kmer_subbed, kmer_subbed.size(), 3);

    nthash.roll();
    nthash.sub({ 3, 4 }, { 'C', 'G' });
    nthash_subbed.roll();
    assert(nthash.get_hash_num() == nthash_subbed.get_hash_num());
    size_t i;
    for (i = 0; i < nthash.get_hash_num(); ++i) {
      assert(nthash.hashes()[i] == nthash_subbed.hashes()[i]);
    }
    assert(i == 3);
  }

  {
    std::cerr << "Testing reverse complement" << std::endl;
    /* Reverse complement of kmer*/
    std::string rc_kmer = "AGACTCAGTCCAGTGTACGT";

    btllib::NtHash nthash(kmer, 20, 3);
    btllib::NtHash nthash_rc(rc_kmer, 20, 3);

    nthash.roll();
    nthash_rc.roll();
    assert(nthash.get_hash_num() == nthash_rc.get_hash_num());
    size_t i;
    for (i = 0; i < nthash.get_hash_num(); ++i) {
      assert(nthash.hashes()[i] == nthash_rc.hashes()[i]);
    }
    assert(i == 3);
  }

  {
    std::cerr << "Testing rolling hash values" << std::endl;
    btllib::NtHash nthash(kmer, 18, 3);

    /* 18-mers of kmer*/
    std::string kmer1 = "ACGTACACTGGACTGAGT";
    std::string kmer2 = "CGTACACTGGACTGAGTC";
    std::string kmer3 = "GTACACTGGACTGAGTCT";

    std::vector<btllib::NtHash> nthash_vector = {
      btllib::NtHash(kmer1, kmer1.size(), nthash.get_hash_num()),
      btllib::NtHash(kmer2, kmer2.size(), nthash.get_hash_num()),
      btllib::NtHash(kmer3, kmer3.size(), nthash.get_hash_num())
    };

    size_t i;
    for (i = 0; nthash.roll() && nthash_vector[i].roll(); ++i) {
      for (size_t j = 0; j < nthash.get_hash_num(); ++j) {
        assert(nthash.hashes()[j] == nthash_vector[i].hashes()[j]);
      }
    }
    assert(i == 3);
  }

  {
    std::cerr << "Testing spaced seeds" << std::endl;
    std::vector<std::string> seeds = { "11111100000000111111",
                                       "11111111000011111111" };

    /* Point Mutations of Kmer*/
    std::string kmerM1 = "ACGTACACTTGACTGAGTCT";
    std::string kmerM2 = "ACGTACACTGTACTGAGTCT";
    std::string kmerM3 = "ACGTACACTGCACTGAGTCT";
    assert(kmerM1.size() == seeds[0].size());
    assert(kmerM1.size() == seeds[1].size());

    btllib::SeedNtHash seed_nthash(kmer, kmer.size(), seeds, 2);

    std::vector<btllib::SeedNtHash> seed_nthash_vector = {
      btllib::SeedNtHash(kmerM1,
                         seed_nthash.get_k(),
                         seeds,
                         seed_nthash.get_hash_num_per_seed()),
      btllib::SeedNtHash(kmerM2,
                         seed_nthash.get_k(),
                         seeds,
                         seed_nthash.get_hash_num_per_seed()),
      btllib::SeedNtHash(
        kmerM3, seed_nthash.get_k(), seeds, seed_nthash.get_hash_num_per_seed())
    };
    assert(seed_nthash.get_hash_num() == seeds.size() * 2);
    assert(seed_nthash.get_hash_num() == seed_nthash_vector[0].get_hash_num());

    seed_nthash.roll();
    size_t i;
    for (i = 0; i < seed_nthash_vector.size() && seed_nthash_vector[i].roll();
         i++) {
      for (size_t j = 0; j < seed_nthash.get_hash_num(); j++) {
        assert(seed_nthash.hashes()[j] == seed_nthash_vector[i].hashes()[j]);
      }
    }
    assert(i == 3);
  }

  {
    std::cerr << "Testing RNA" << std::endl;
    btllib::NtHash dna_nthash(kmer, 20, 3);

    std::string rna_kmer = "ACGUACACUGGACUGAGUCU";
    btllib::NtHash rna_nthash(kmer, 20, 3);

    dna_nthash.roll();
    rna_nthash.roll();
    size_t i;
    for (i = 0; i < dna_nthash.get_hash_num(); ++i) {
      assert(dna_nthash.hashes()[i] == rna_nthash.hashes()[i]);
    }
    assert(i == 3);
  }

  return 0;
}