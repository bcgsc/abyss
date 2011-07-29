/*
 *  Copyright (c) 2010 Yasuo Tabei
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

#ifndef _FMINDEX_HPP_
#define _FMINDEX_HPP_

#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <set>
#include <cmath>
#include "stdint.h"
#include <ctime>
#include <limits.h>
#include "wat_array.h"

#define UINT8_MAX (0x100 - 1)

class FMIndex {
  int                              percent;
  uint8_t                          alphaSize;
  std::vector<uint32_t>            cf;
  std::map<unsigned char, uint8_t> mapping;
  std::map<uint8_t, unsigned char> rmapping;
  std::vector<uint32_t>            sampledSA;
  std::vector<std::vector<int> >   cols;
  wat_array::WatArray              wa;

  void calculateStatistics(const std::vector<uint8_t> &s);
  int  buildBWT(const std::vector<uint8_t> &s, const std::vector<uint32_t> &sa, std::vector<uint64_t> &bwt);
  int  buildSA(const std::vector<uint8_t> &s, std::vector<uint32_t> &sa);
  int  buildSampledSA(const std::vector<uint8_t> &s, const std::vector<uint32_t> &sa);
  int  buildWaveletTree(const std::vector<uint64_t> &bwt);
  void searchHamming(const std::vector<uint8_t> &qs, int i, size_t sp, size_t ep, int d, int dist, std::vector<std::pair<uint64_t, uint64_t> > &res);
  void searchEdit(const std::vector<uint8_t> &qs, int count, int sp, int ep, int dist, std::vector<std::pair<uint64_t, uint64_t> > &res);
public:
  int    load(std::istream& is);
  int    save(std::ostream& os);
  int    read(const char *fname, std::vector<uint8_t> &s);
  int    readQstring(const char *fname, std::vector<std::vector<uint8_t> > &qs);
  int    buildFmIndex(const char *fnmae, int _percent);
  void   search(const std::vector<uint8_t> &qs, std::pair<uint64_t, uint64_t> &res);
  void   searchHamming(const std::vector<uint8_t> &qs, int dist, std::vector<std::pair<uint64_t, uint64_t> > &res);
  void   searchEdit(const std::vector<uint8_t> &qs, int dist, std::vector<std::pair<uint64_t, uint64_t> > &res);
  size_t locate(uint64_t i);
};

#endif // _FMINDEX_H_
