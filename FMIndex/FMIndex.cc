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
#include "sais.h"
#include <climits>
#include <fstream>

using namespace std;

int FMIndex::read(const char *fname, vector<uint8_t> &s) {
  ifstream ifs(fname);
  string line;
  while (getline(ifs, line)) {
    for (size_t i = 0; i < line.size(); i++) {
      if (mapping.find(line[i]) == mapping.end()) {
	mapping[line[i]]    = alphaSize;
	rmapping[alphaSize] = line[i];
	if (++alphaSize == 0)
	  cerr << "warning: the variety of characters exceeds 255." << endl;
      }
      s.push_back(mapping[line[i]]);
    }
    s.push_back(mapping['\n']);
  }
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
  size_t divSeq = seqLen*((float)percent/100.f);
  divSeq = seqLen/divSeq;

  for (size_t i = 0; i < seqLen; i++) {
    if (i%divSeq == 0)
      sampledSA.push_back(sa[i]);
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
  cf.resize(UINT8_MAX + 1);
  cf[0] = 0;
  for (size_t i = 1; i <= UINT8_MAX; i++) {
    cf[i]     = tmpCf[i-1];
    tmpCf[i] += tmpCf[i-1];
  }
}

size_t FMIndex::locate(uint64_t i) const {
  size_t bsize = wa.length();
  size_t divSeq = bsize*((float)percent/100.f);
  divSeq = bsize/divSeq;

  size_t j = i;
  size_t t = 0;
  while (j % divSeq != 0) {
    j = cf[wa.Lookup(j)] + wa.Rank(wa.Lookup(j), j+1) - 1;
    t++;
  }
  if (sampledSA[j/divSeq] + t >= bsize) {
    return (size_t)sampledSA[j/divSeq] + t - bsize;
  }
  return (size_t)sampledSA[j/divSeq] + t;
}

int FMIndex::save(ostream& os)  {
  os.write((const char*)(&percent), sizeof(percent));
  os.write((const char*)(&alphaSize), sizeof(alphaSize));
  size_t cfSize = cf.size();
  os.write((const char*)(&cfSize), sizeof(cfSize));
  size_t sampledSASize = sampledSA.size();
  os.write((const char*)(&sampledSASize), sizeof(sampledSASize));
  os.write((const char*)(&cf[0]), sizeof(cf[0])*cf.size());
  os.write((const char*)(&sampledSA[0]), sizeof(sampledSA[0])*sampledSA.size());

  uint32_t count = 0;
  vector<unsigned char> key;
  vector<uint8_t> val;
  for (std::map<unsigned char, uint8_t>::iterator it = mapping.begin(); it != mapping.end(); it++) {
    key.push_back(it->first);
    val.push_back(it->second);
    count++;
  }
  os.write((const char*)(&count), sizeof(uint32_t));
  os.write((const char*)(&key[0]), sizeof(key[0])*count);
  os.write((const char*)(&val[0]), sizeof(val[0])*count);

  wa.Save(os);

  return 0;
}

int FMIndex::load(istream& is) {
  is.read((char*)(&percent), sizeof(percent));
  is.read((char*)(&alphaSize), sizeof(alphaSize));
  size_t cfSize;
  is.read((char*)(&cfSize), sizeof(cfSize));
  size_t sampledSASize;
  is.read((char*)(&sampledSASize), sizeof(sampledSASize));
  cf.resize(cfSize);
  is.read((char*)(&cf[0]), sizeof(cf[0])*cfSize);
  sampledSA.resize(sampledSASize);
  is.read((char*)(&sampledSA[0]), sizeof(sampledSA[0])*sampledSASize);

  uint32_t              count;
  vector<unsigned char> key;
  vector<uint8_t>       val;
  is.read((char*)(&count), sizeof(count));
  key.resize(count);
  is.read((char*)(&key[0]), sizeof(key[0])*count);
  val.resize(count);
  is.read((char*)(&val[0]), sizeof(val[0])*count);
  mapping.clear();
  rmapping.clear();
  for (uint32_t i = 0; i < count; i++) {
    mapping[key[i]]  = val[i];
    rmapping[val[i]] = key[i];
  }
  wa.Load(is);

  return 0;
}

int FMIndex::buildWaveletTree(const vector<uint64_t> &bwt) {
  wa.Init(bwt);

  return 0;
}

int FMIndex::buildFmIndex(const char *fname, int _percent) {
  percent = _percent;

  cerr << "start reading the input-file" << endl;
  vector<uint8_t> s;
  read(fname, s);

  cerr << "alphabet size:" << (int)alphaSize << endl;

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
