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

#include "FMIndex.h"
#include "IOUtil.h"
#include "sais.h"
#include <algorithm>
#include <cassert>
#include <climits>
#include <fstream>

using namespace std;

int FMIndex::read(const char *fname, vector<uint8_t> &s)
{
	// Read the file.
	ifstream in(fname);
	assert_good(in, fname);
	in.seekg(0, ios::end);
	s.resize(in.tellg());
	in.seekg(0, ios::beg);
	assert_good(in, fname);
	in.read(reinterpret_cast<char*>(s.data()), s.size());
	assert_good(in, fname);
	assert((size_t)in.gcount() == s.size());

	// Translate the alphabet.
	if (m_alphaSize == 0)
		setAlphabet(s.begin(), s.end());
	transform(s.begin(), s.end(), s.begin(), Translate(*this));
	replace(s.begin(), s.end(), UCHAR_MAX, 1);

	// Add the sentinel character.
	s.push_back(0);

	return 0;
}

int FMIndex::buildSA(const vector<uint8_t> &s, vector<uint32_t> &sa) {
  sa.resize(s.size());
  if (saisxx(s.begin(), sa.begin(), (int)s.size(), 0x100) == -1) {
    cerr << "suffix array construction error" << endl;
    return 1;
  }

  return 0;
}

int FMIndex::buildBWT(const vector<uint8_t> &s, const vector<uint32_t> &sa, vector<uint64_t> &bwt) {
  size_t seqLen = s.size();
  bwt.resize(seqLen);
  for (size_t i = 0; i < seqLen; i++) {
    if (sa[i] > 0) bwt[i] = s[sa[i] - 1];
    else           bwt[i] = s[seqLen - 1];
  }

  return 0;
}

int FMIndex::buildSampledSA(const vector<uint8_t> &s, const vector<uint32_t> &sa) {
  size_t seqLen = s.size();
  size_t divSeq = seqLen*((float)m_percent/100.f);
  divSeq = seqLen/divSeq;

  for (size_t i = 0; i < seqLen; i++) {
    if (i%divSeq == 0)
      m_sampledSA.push_back(sa[i]);
  }
  return 0;
}

void FMIndex::calculateStatistics(const vector<uint8_t> &s) {
  const unsigned UINT8_MAX = 255;
  size_t seqLen = s.size();
  vector<uint32_t> tmpCf(UINT8_MAX + 1);
  for (size_t i = 0; i < seqLen; i++) {
    tmpCf[s[i]]++;
  }
  m_cf.resize(UINT8_MAX + 1);
  m_cf[0] = 0;
  for (size_t i = 1; i <= UINT8_MAX; i++) {
    m_cf[i] = tmpCf[i-1];
    tmpCf[i] += tmpCf[i-1];
  }
}

size_t FMIndex::locate(uint64_t i) const {
  size_t bsize = m_wa.length();
  size_t divSeq = bsize*((float)m_percent/100.f);
  divSeq = bsize/divSeq;

  size_t j = i;
  size_t t = 0;
  while (j % divSeq != 0) {
    unsigned c = m_wa.Lookup(j);
    j = m_cf[c] + m_wa.Rank(c, j+1) - 1;
    t++;
  }
  if (m_sampledSA[j/divSeq] + t >= bsize) {
    return (size_t)m_sampledSA[j/divSeq] + t - bsize;
  }
  return (size_t)m_sampledSA[j/divSeq] + t;
}

int FMIndex::save(ostream& os)  {
  os.write((const char*)(&m_percent), sizeof(m_percent));
  os.write((const char*)(&m_alphaSize), sizeof(m_alphaSize));
  size_t cfSize = m_cf.size();
  os.write((const char*)(&cfSize), sizeof(cfSize));
  size_t sampledSASize = m_sampledSA.size();
  os.write((const char*)(&sampledSASize), sizeof(sampledSASize));
  os.write((const char*)(&m_cf[0]), sizeof(m_cf[0])*m_cf.size());
  os.write((const char*)&m_sampledSA[0],
		  sizeof m_sampledSA[0] * m_sampledSA.size());

  vector<unsigned char> key;
  vector<uint8_t> val;
  for (unsigned i = 0; i < m_mapping.size(); ++i) {
    if (m_mapping[i] != UCHAR_MAX) {
      key.push_back(i);
      val.push_back(m_mapping[i]);
    }
  }
  uint32_t count = key.size();
  os.write((const char*)(&count), sizeof(uint32_t));
  os.write((const char*)(&key[0]), sizeof(key[0])*count);
  os.write((const char*)(&val[0]), sizeof(val[0])*count);

  m_wa.Save(os);

  return 0;
}

int FMIndex::load(istream& is) {
  is.read((char*)(&m_percent), sizeof(m_percent));
  is.read((char*)(&m_alphaSize), sizeof(m_alphaSize));
  size_t cfSize;
  is.read((char*)(&cfSize), sizeof(cfSize));
  size_t sampledSASize;
  is.read((char*)(&sampledSASize), sizeof(sampledSASize));
  m_cf.resize(cfSize);
  is.read((char*)(&m_cf[0]), sizeof(m_cf[0])*cfSize);
  m_sampledSA.resize(sampledSASize);
  is.read((char*)&m_sampledSA[0],
		  sizeof m_sampledSA[0] *sampledSASize);

  uint32_t              count;
  vector<unsigned char> key;
  vector<uint8_t>       val;
  is.read((char*)(&count), sizeof(count));
  key.resize(count);
  is.read((char*)(&key[0]), sizeof(key[0])*count);
  val.resize(count);
  is.read((char*)(&val[0]), sizeof(val[0])*count);
  m_mapping.clear();
  for (uint32_t i = 0; i < count; i++) {
    unsigned k = key[i];
    unsigned v = val[i];
    if (k >= m_mapping.size())
      m_mapping.resize(k + 1, UCHAR_MAX);
    assert(m_mapping[k] == UCHAR_MAX);
    m_mapping[k] = v;
  }
  m_wa.Load(is);

  return 0;
}

int FMIndex::buildWaveletTree(const vector<uint64_t> &bwt) {
  m_wa.Init(bwt);

  return 0;
}

int FMIndex::buildFmIndex(const char *fname, int percent) {
  m_percent = percent;

  cerr << "start reading the input-file" << endl;
  vector<uint8_t> s;
  read(fname, s);

  cerr << "alphabet size:" << (int)m_alphaSize << endl;

  double sTime = clock();
  vector<uint32_t> sa;
  cerr << "build SA" << endl;
  buildSA(s, sa);

  cerr << "calculate statistics" << endl;
  calculateStatistics(s);

  vector<uint64_t> bwt;
  cerr << "build BWT" << endl;
  buildBWT(s, sa, bwt);

  cerr << "build WaveletTree" << endl;
  buildWaveletTree(bwt);

  cerr << "build sampledSA" << endl;
  buildSampledSA(s, sa);

  double eTime = clock();

  fprintf(stderr, "cpu time:%f\n", (eTime - sTime)/CLOCKS_PER_SEC);

  return 0;
}
