#include "btllib/indexlr.hpp"
#include "btllib/bloom_filter.hpp"

#include <fstream>
#include <sstream>

int
main()
{
  btllib::Indexlr indexlr("../tests/indexlr.fa", 100, 5, 0);
  btllib::Indexlr indexlr2("../tests/indexlr.fq",
                           100,
                           5,
                           btllib::Indexlr::Flag::ID |
                             btllib::Indexlr::Flag::BX |
                             btllib::Indexlr::Flag::SEQ);

  std::ifstream correct_output_file("../tests/indexlr.fa.correct");
  std::string correct_output;
  correct_output_file.seekg(0, std::ios::end);
  correct_output.reserve(correct_output_file.tellg());
  correct_output_file.seekg(0, std::ios::beg);

  correct_output.assign(std::istreambuf_iterator<char>(correct_output_file),
                        std::istreambuf_iterator<char>());

  std::ifstream correct_output_file2("../tests/indexlr.fq.correct");
  std::string correct_output2;
  correct_output_file2.seekg(0, std::ios::end);
  correct_output2.reserve(correct_output_file2.tellg());
  correct_output_file2.seekg(0, std::ios::beg);

  correct_output2.assign(std::istreambuf_iterator<char>(correct_output_file2),
                         std::istreambuf_iterator<char>());

  std::stringstream ss;
  std::stringstream ss2;

  std::cerr << "Testing without Bloom filters" << std::endl;
  decltype(indexlr)::Record record;
  bool success_indexlr = false, success_indexlr2 = false;
  for (int i = 0;; i++) {
    if ((success_indexlr = (record = indexlr.get_minimizers()))) {
      if (i > 0) {
        ss << '\n';
      }
      ss << record.id << '\t';
      int j = 0;
      for (const auto& min : record.minimizers) {
        if (j > 0) {
          ss << ' ';
        }
        ss << min.out_hash;
        j++;
      }
    }
    if ((success_indexlr2 = (record = indexlr2.get_minimizers()))) {
      if (i > 0) {
        ss2 << '\n';
      }
      ss2 << record.id << '\t' << record.barcode << '\t';
      int j = 0;
      for (const auto& min : record.minimizers) {
        if (j > 0) {
          ss2 << ' ';
        }
        ss2 << min.out_hash << ':' << min.pos << ':'
            << (min.forward ? '+' : '-') << ':' << min.seq;
        j++;
      }
    }
    if (!success_indexlr && !success_indexlr2) {
      break;
    }
  }

  if (ss.str() != correct_output) {
    std::cerr << "Correct:\n\n"
              << correct_output << "\n\nActual:\n\n"
              << ss.str() << "\n\n";
  }
  assert(ss.str() == correct_output);

  if (ss2.str() != correct_output2) {
    std::cerr << "Correct:\n\n"
              << correct_output2 << "\n\nActual:\n\n"
              << ss2.str() << "\n\n";
  }
  assert(ss2.str() == correct_output2);

  std::cerr << "Testing with Bloom filters" << std::endl;
  btllib::BloomFilter filter_in_bf(1024 * 1024 * 32, 1);
  btllib::BloomFilter filter_out_bf(1024 * 1024 * 32, 1);

  std::vector<uint64_t> filter_in_hashes = { 430447521414431149ULL,
                                             3146270839399521840ULL,
                                             161808173335193798ULL };
  std::vector<uint64_t> filter_out_hashes = { 1672947938795563804ULL,
                                              2858314356342515870ULL,
                                              1712341822067595113ULL };

  for (const auto h : filter_in_hashes) {
    filter_in_bf.insert({ h });
  }
  for (const auto h : filter_out_hashes) {
    filter_out_bf.insert({ h });
  }

  btllib::Indexlr indexlr3("../tests/indexlr.fq",
                           100,
                           5,
                           btllib::Indexlr::Flag::FILTER_IN,
                           3,
                           true,
                           filter_in_bf);
  size_t mins_found = 0;
  while ((record = indexlr3.get_minimizers())) {
    for (const auto& min : record.minimizers) {
      bool found = false;
      for (const auto h : filter_in_hashes) {
        if (min.min_hash == h) {
          found = true;
          break;
        }
      }
      assert(found);
      mins_found++;
    }
  }
  assert(mins_found >= filter_in_hashes.size());

  btllib::Indexlr indexlr4("../tests/indexlr.fq",
                           100,
                           5,
                           btllib::Indexlr::Flag::FILTER_OUT,
                           3,
                           true,
                           filter_out_bf);
  mins_found = 0;
  while ((record = indexlr4.get_minimizers())) {
    for (const auto& min : record.minimizers) {
      for (const auto h : filter_out_hashes) {
        assert(min.min_hash != h);
      }
      mins_found++;
    }
  }
  assert(mins_found >= filter_in_hashes.size());

  btllib::Indexlr indexlr5("../tests/indexlr.fq",
                           100,
                           5,
                           btllib::Indexlr::Flag::FILTER_IN |
                             btllib::Indexlr::Flag::FILTER_OUT,
                           3,
                           true,
                           filter_in_bf,
                           filter_out_bf);
  mins_found = 0;
  while ((record = indexlr5.get_minimizers())) {
    for (const auto& min : record.minimizers) {
      bool found = false;
      for (const auto h : filter_in_hashes) {
        if (min.min_hash == h) {
          found = true;
          break;
        }
      }
      assert(found);
      for (const auto h : filter_out_hashes) {
        assert(min.min_hash != h);
      }
      mins_found++;
    }
  }
  assert(mins_found >= filter_in_hashes.size());

  return 0;
}