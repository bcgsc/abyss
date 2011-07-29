/*
 *  Copyright (c) 2010 Daisuke Okanohara
 *
 *   Redistribution and use in source and binary forms, with or without
 *   modification, are permitted provided that the following conditions
 *   are met:
 *
 *   1. Redistributions of source code must retain the above Copyright
 *      notice, this list of conditions and the following disclaimer.
 *
 *   2. Redistributions in binary form must reproduce the above Copyright
 *      notice, this list of conditions and the following disclaimer in the
 *      documentation and/or other materials provided with the distribution.
 *
 *   3. Neither the name of the authors nor the names of its contributors
 *      may be used to endorse or promote products derived from this
 *      software without specific prior written permission.
 */

#ifndef WATARRAY_WATARRAY_HPP_
#define WATARRAY_WATARRAY_HPP_

#include <vector>
#include <queue>
#include <stdint.h>
#include <iostream>
#include <cassert>
#include "bit_array.h"

namespace wat_array {

/**
 Wavelet Tree Array Libarary for array processing.

 Input: A[0...n], 0 <= A[i] < k,
 Space: n log_2 k bits

 Support many queries in O(log k) time and constant for n.
 */

struct ListResult{
  ListResult(uint64_t c, uint64_t freq) : c(c), freq(freq){}
  uint64_t c;
  uint64_t freq;
  bool operator < (const ListResult& lr) const {
    if (c != lr.c) return c < lr.c;
    return freq < lr.freq;
  }
};


class WatArray{
public:
  /**
   * Constructor
   */
  WatArray();

  /**
   * Destructor
   */
  ~WatArray();

  /**
   * Initialize an index from an array
   * @param An array to be initialized
   */
  void Init(const std::vector<uint64_t>& array);

  /**
   * Clear and release the resouces
   */
  void Clear();

  /**
   * Lookup A[pos]
   * @param pos the position
   * @return return A[pos] if found, or return NOT_FOUND if pos >= length
   */
  uint64_t Lookup(uint64_t pos) const;

  /**
   * Compute the rank = the frequency of a character 'c' in the prefix of the array A[0...pos)
   * @param c Character to be examined
   * @param pos The position of the prefix (not inclusive)
   * @return The frequency of a character 'c' in the prefix of the array A[0...pos)
   *         or NOT_FOUND if c >= alphabet_num or pos > length
   */
  uint64_t Rank(uint64_t c, uint64_t pos) const;

  /**
   * Compute the select = the position of the (rank+1)-th occurence of 'c' in the array.
   * @param c Character to be examined
   * @param rank The rank of the character
   * @return The position of the (rank+1)-th occurence of 'c' in the array.
   *         or NOT_FOUND if c >= alphabet_num or rank > Freq(c)
   */
  uint64_t Select(uint64_t c, uint64_t rank) const;

  /**
   * Compute the frequency of characters c' < c in the subarray A[0...pos)
   * @param c The upper bound of the characters
   * @param pos The position of the end of the prefix (not inclusive)
   * @return The frequency of characters c' < c in the prefix of the array A[0...pos)
             or NOTFOUND if c > alphabet_num or pos > length
   */
  uint64_t RankLessThan(uint64_t c, uint64_t pos) const;

  /**
   * Compute the frequency of characters c' > c in the subarray A[0...pos)
   * @param c The lower bound of the characters
   * @param pos The position of the end of the prefix (not inclusive)
   * @return The frequency of characters c' < c in the prefix of the array A[0...pos)
             or NOTFOUND if c > alphabet_num or pos > length
   */
  uint64_t RankMoreThan(uint64_t c, uint64_t pos) const;

  /**
   * Compute the frequency of characters c' < c, c'=c, and c' > c, in the subarray A[0...pos)
   * @param c The character
   * @param pos The position of the end of the prefix (not inclusive)
   * @param rank The frefquency of c in A[0...pos)
   * @param rank_less_than The frequency of c' < c in A[0...pos)
   * @param rank_more_than The frequency of c' > c in A[0...pos)
    */
  void RankAll(uint64_t c, uint64_t pos, uint64_t& rank,
	       uint64_t& rank_less_than, uint64_t& rank_more_than) const;

  /**
   * Compute the frequency of characters min_c <= c' < max_c in the subarray A[beg_pos ... end_pos)
   * @param min_c The smallerest character to be examined
   * @param max_c The uppker bound of the character to be examined
   * @param beg_pos The beginning position of the array (inclusive)
   * @param end_pos The ending position of the array (not inclusive)
   * @return The frequency of characters min_c <= c < max_c in the subarray A[beg_pos .. end_pos)
             or NOTFOUND if max_c > alphabet_num or end_pos > length
   */
  uint64_t FreqRange(uint64_t min_c, uint64_t max_c, uint64_t beg_pos, uint64_t end_pos) const;

  /**
   * Range Max Query
   * @param beg_pos The beginning position
   * @param end_pos The ending position
   * @param pos The position where the largest value appeared in the subarray A[beg_pos .. end_pos)
                If there are many items having the largest values, the smallest pos will be reporeted
   * @param val The largest value appeared in the subarray A[beg_pos ... end_pos)
   */
  void MaxRange(uint64_t beg_pos, uint64_t end_pos, uint64_t& pos, uint64_t& val) const;

  /**
   * Range Min Query
   * @param beg_pos The beginning position
   * @param end_pos The ending position
   * @param pos The position where the smallest value appeared in the subarray A[beg_pos .. end_pos)
                If there are many items having the smalles values, the smallest pos will be reporeted
   * @param val The smallest value appeared in the subarray A[beg_pos ... end_pos)
   */
  void MinRange(uint64_t beg_pos, uint64_t end_pos, uint64_t& pos, uint64_t& val) const;

  /**
   * Range Quantile Query, Return the K-th smallest value in the subarray
   * @param beg_pos The beginning position
   * @param end_pos The ending position
   * @param k The order (should be smaller than end_pos - beg_pos).
   * @param pos The position where the k-th largest value appeared in the subarray A[beg_pos .. end_pos)
                If there are many items having the k-th largest values, the smallest pos will be reporeted
   * @param val The k-th largest value appeared in the subarray A[beg_pos ... end_pos)
   */
  void QuantileRange(uint64_t beg_pos, uint64_t end_pos, uint64_t k, uint64_t& pos, uint64_t& val) const;

  /**
   * List the distinct characters appeared in A[beg_pos ... end_pos) from most frequent ones
   */
  void ListModeRange(uint64_t min_c, uint64_t max_c, uint64_t beg_pos, uint64_t end_pos, uint64_t num, std::vector<ListResult>& res) const;

  /**
   * List the distinct characters in A[beg_pos ... end_pos) min_c <= c < max_c  from smallest ones
   * @param min_c The smallerest character to be examined
   * @param max_c The uppker bound of the character to be examined
   * @param beg_pos The beginning position of the array (inclusive)
   * @param end_pos The ending positin of the array (not inclusive)
   * @param num The maximum number of reporting results.
   * @param res The distinct chracters in the A[beg_pos ... end_pos) from smallest ones.
   *            Each item consists of c:character and freq: frequency of c.
   */
  void ListMinRange(uint64_t min_c, uint64_t max_c, uint64_t beg_pos, uint64_t end_pos, uint64_t num, std::vector<ListResult>& res) const;

  /**
   * List the distinct characters appeared in A[beg_pos ... end_pos) from largest ones
   * @param min_c The smallerest character to be examined
   * @param max_c The uppker bound of the character to be examined
   * @param beg_pos The beginning position of the array (inclusive)
   * @param end_pos The ending positin of the array (not inclusive)
   * @param num The maximum number of reporting results.
   * @param res The distinct chracters in the A[beg_pos ... end_pos) from largestx ones.
   *            Each item consists of c:character and freq: frequency of c.
   */
  void ListMaxRange(uint64_t min_c, uint64_t max_c, uint64_t beg_pos, uint64_t end_pos, uint64_t num, std::vector<ListResult>& res) const;

  /**
   * Compute the frequency of the character c
   * @param c The character to be examined
   * param Return the frequency of c in the array.
   */
  uint64_t Freq(uint64_t c) const;

  /**
   * Compute the frequency of the characters
   * @param min_c The minimum character
   * @param max_c The maximum character
   * param Return the frequency of min_c <= c < max_c in the array.
   */
  uint64_t FreqSum(uint64_t min_c, uint64_t max_c) const;

  /**
   * Return the number of alphabets in the array
   * @return The number of alphabet in the array
   */
  uint64_t alphabet_num() const;

  /**
   * Return the length of the array
   * @return The length of the array
   */
  uint64_t length() const;

  /**
   * Save the current status to a stream
   * @param os The output stream where the data is saved
   */
  void Save(std::ostream& os) const;

  /**
   * Load the current status from a stream
   * @param is The input stream where the status is saved
   */
  void Load(std::istream& is);

private:
  uint64_t GetAlphabetNum(const std::vector<uint64_t>& array) const;
  uint64_t Log2(uint64_t x) const;
  uint64_t PrefixCode(uint64_t x, uint64_t len, uint64_t total_len) const;
  static uint64_t GetMSB(uint64_t x, uint64_t pos, uint64_t len);
  static uint64_t GetLSB(uint64_t x, uint64_t pos);
  void SetArray(const std::vector<uint64_t>& array);
  void SetOccs(const std::vector<uint64_t>& array);
  void GetBegPoses(const std::vector<uint64_t>& array, uint64_t alpha_bit_num,
		   std::vector<std::vector<uint64_t> >& beg_poses) const;

  struct QueryOnNode{
    QueryOnNode(uint64_t beg_node, uint64_t end_node, uint64_t beg_pos, uint64_t end_pos,
		uint64_t depth, uint64_t prefix_char) :
      beg_node(beg_node), end_node(end_node), beg_pos(beg_pos), end_pos(end_pos),
      depth(depth), prefix_char(prefix_char) {}
    uint64_t beg_node;
    uint64_t end_node;
    uint64_t beg_pos;
    uint64_t end_pos;
    uint64_t depth;
    uint64_t prefix_char;
    void print() {
      std::cout << beg_node << " " << end_node << " "
		<< beg_pos  << " " << end_pos  << " "
		<< depth << " " ;
      for (uint64_t i = 0; i < depth; ++i){
	std::cout << GetMSB(prefix_char, i, depth);
      }
      std::cout << std::endl;
    }
  };

  class ListModeComparator;
  class ListMinComparator;
  class ListMaxComparator;

  template <class Comparator>
  void ListRange(uint64_t min_c,   uint64_t max_c,
		 uint64_t beg_pos, uint64_t end_pos,
		 uint64_t num, std::vector<ListResult>& res) const {
    res.clear();
    if (end_pos > length_ || beg_pos >= end_pos) return;

    std::priority_queue<QueryOnNode, std::vector<QueryOnNode>, Comparator> qons;
    qons.push(QueryOnNode(0, length_, beg_pos, end_pos, 0, 0));

    while (res.size() < num && !qons.empty()){
      QueryOnNode qon = qons.top();
      qons.pop();
      if (qon.depth >= alphabet_bit_num_){
	res.push_back(ListResult(qon.prefix_char, qon.end_pos - qon.beg_pos));
      } else {
	std::vector<QueryOnNode> next;
	ExpandNode(min_c, max_c, qon, next);
	for (size_t i = 0; i < next.size(); ++i){
	  qons.push(next[i]);
	}
      }
    }
  }

  bool CheckPrefix(uint64_t prefix, uint64_t depth, uint64_t min_c, uint64_t max_c) const;
  void ExpandNode(uint64_t min_c, uint64_t max_c,
		  const QueryOnNode& qon, std::vector<QueryOnNode>& next) const;

  std::vector<BitArray> bit_arrays_;
  BitArray occs_;

  uint64_t alphabet_num_;
  uint64_t alphabet_bit_num_;
  uint64_t length_;
};



}


#endif // WASEQ_WASEQ_HPP_
