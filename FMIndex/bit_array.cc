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

#include "bit_array.h"
#include "BitUtil.h" // for popcount
#include <cassert>

namespace wat_array {

BitArray::BitArray() : length_(0), one_num_(0){
}

BitArray::BitArray(uint64_t length){
  Init(length);
}

BitArray::~BitArray(){
}

uint64_t BitArray::length() const {
  return length_;
}

uint64_t BitArray::one_num() const{
  return one_num_;
}

void BitArray::Init(uint64_t length){
  length_    = length;
  one_num_ = 0;
  uint64_t block_num = (length + BLOCK_BITNUM - 1) / BLOCK_BITNUM;
  bit_blocks_.resize(block_num);
}

void BitArray::Clear(){
  std::vector<uint64_t>().swap(bit_blocks_);
  std::vector<uint64_t>().swap(rank_tables_);
  length_ = 0;
  one_num_ = 0;
}

void BitArray::Build() {
  one_num_ = 0;
  uint64_t table_num = ((bit_blocks_.size() + TABLE_INTERVAL - 1) / TABLE_INTERVAL) + 1;
  rank_tables_.resize(table_num);
  for (size_t i = 0; i < bit_blocks_.size(); ++i){
    if ((i % TABLE_INTERVAL) == 0){
      rank_tables_[i/TABLE_INTERVAL] = one_num_;
    }
    one_num_ += PopCount(bit_blocks_[i]);
  }
  rank_tables_.back() = one_num_;
}

void BitArray::SetBit(uint64_t bit, uint64_t pos) {
  if (!bit) return;
  bit_blocks_[pos / BLOCK_BITNUM] |= (1LLU << (pos % BLOCK_BITNUM));
}

uint64_t BitArray::Rank(uint64_t bit, uint64_t pos) const {
  if (pos > length_) return NOTFOUND;
  if (bit) return RankOne(pos);
  else return pos - RankOne(pos);
}

uint64_t BitArray::Select(uint64_t bit, uint64_t rank) const {
  if (bit){
    if (rank > one_num_) return NOTFOUND;
  } else {
    if (rank > length_ - one_num_) return NOTFOUND;
  }

  uint64_t block_pos = SelectOutBlock(bit, rank);
  uint64_t block = (bit) ? bit_blocks_[block_pos] : ~bit_blocks_[block_pos];
  return block_pos * BLOCK_BITNUM + SelectInBlock(block, rank);
}

uint64_t BitArray::SelectOutBlock(uint64_t bit, uint64_t& rank) const {
  // binary search over tables
  uint64_t left = 0;
  uint64_t right = rank_tables_.size();
  while (left < right){
    uint64_t mid = (left + right) / 2;
    uint64_t length = BLOCK_BITNUM * TABLE_INTERVAL * mid;
    if (GetBitNum(rank_tables_[mid], length, bit) < rank) {
      left = mid+1;
    } else {
      right = mid;
    }
  }

  uint64_t table_ind   = (left != 0) ? left - 1: 0;
  uint64_t block_pos   = table_ind * TABLE_INTERVAL;
  rank -= GetBitNum(rank_tables_[table_ind],
		    block_pos * BLOCK_BITNUM,
		    bit);

  // sequential search over blocks
  for ( ; block_pos < bit_blocks_.size(); ++block_pos){
    uint64_t rank_next= GetBitNum(PopCount(bit_blocks_[block_pos]), BLOCK_BITNUM, bit);
    if (rank <= rank_next){
      break;
    }
    rank -= rank_next;
  }
  return block_pos;
}

uint64_t BitArray::SelectInBlock(uint64_t x, uint64_t rank) {
  uint64_t x1 = x - ((x & 0xAAAAAAAAAAAAAAAALLU) >> 1);
  uint64_t x2 = (x1 & 0x3333333333333333LLU) + ((x1 >> 2) & 0x3333333333333333LLU);
  uint64_t x3 = (x2 + (x2 >> 4)) & 0x0F0F0F0F0F0F0F0FLLU;

  uint64_t pos = 0;
  for (;;  pos += 8){
    uint64_t rank_next = (x3 >> pos) & 0xFFLLU;
    if (rank <= rank_next) break;
    rank -= rank_next;
  }

  uint64_t v2 = (x2 >> pos) & 0xFLLU;
  if (rank > v2) {
    rank -= v2;
    pos += 4;
  }

  uint64_t v1 = (x1 >> pos) & 0x3LLU;
  if (rank > v1){
    rank -= v1;
    pos += 2;
  }

  uint64_t v0  = (x >> pos) & 0x1LLU;
  if (v0 < rank){
    rank -= v0;
    pos += 1;
  }

  return pos;
}

uint64_t BitArray::Lookup(uint64_t pos) const {
  return (bit_blocks_[pos / BLOCK_BITNUM] >> (pos % BLOCK_BITNUM)) & 1LLU;
}


uint64_t BitArray::RankOne(uint64_t pos) const {
  uint64_t block_ind = pos / BLOCK_BITNUM;
  uint64_t table_ind = block_ind / TABLE_INTERVAL;
  assert(table_ind < rank_tables_.size());

  uint64_t rank = rank_tables_[table_ind];
  for (uint64_t i = table_ind * TABLE_INTERVAL; i < block_ind; ++i){
    rank += PopCount(bit_blocks_[i]);
  }
  rank += PopCountMask(bit_blocks_[block_ind], pos % BLOCK_BITNUM);
  return rank;
}

/** Return the Hamming weight of x. */
uint64_t BitArray::PopCount(uint64_t x)
{
	return popcount(x);
}

uint64_t BitArray::PopCountMask(uint64_t x, uint64_t offset) {
  if (offset == 0) return 0;
  return PopCount(x & ((1LLU << offset) - 1));
}

uint64_t BitArray::GetBitNum(uint64_t one_num, uint64_t num, uint64_t bit) {
  if (bit) return one_num;
  else return num - one_num;
}

void BitArray::PrintForDebug(std::ostream& os) const {
  for (uint64_t i = 0; i < length_;  ++i){
    if (Lookup(i)) os << "1";
    else           os << "0";
    if (((i+1) % 8) == 0) {
      os << " ";
    }
  }
}

void BitArray::Save(std::ostream& os) const{
  os.write((const char*)(&length_), sizeof(length_));
  os.write((const char*)(&bit_blocks_[0]), sizeof(bit_blocks_[0]) * bit_blocks_.size());
}

void BitArray::Load(std::istream& is){
  Clear();
  is.read((char*)(&length_), sizeof(length_));
  Init(length_);
  is.read((char*)(&bit_blocks_[0]), sizeof(bit_blocks_[0]) * bit_blocks_.size());
  Build();
}


}
