/*
 * Copyright (c) 2010 Yasuo Tabei
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above Copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above
 *    Copyright notice, this list of conditions and the following
 *    disclaimer in the documentation and/or other materials provided
 *    with the distribution.
 *
 * 3. Neither the name of the authors nor the names of its
 *    contributors may be used to endorse or promote products derived
 *    from this software without specific prior written permission.
 */

#include "FMIndex.h"
#include "IOUtil.h"
#include "sais.h"
#include <algorithm>
#include <cassert>
#include <climits>
#include <cstdlib>
#include <ctime> // for clock
#include <fstream>
#include <iostream>
#include <stdint.h>
#include <vector>

using namespace std;

/** Read the sequence to be indexed. */
void FMIndex::read(const char* path, vector<uint8_t>& s)
{
	// Read the file.
	ifstream in(path);
	assert_good(in, path);
	in.seekg(0, ios::end);
	s.resize(in.tellg());
	in.seekg(0, ios::beg);
	assert_good(in, path);
	in.read((char*)s.data(), s.size());
	assert_good(in, path);
	assert((size_t)in.gcount() == s.size());

	// Translate the alphabet.
	if (m_alphaSize == 0)
		setAlphabet(s.begin(), s.end());
	transform(s.begin(), s.end(), s.begin(), Translate(*this));
	replace(s.begin(), s.end(), UCHAR_MAX, 1);

	// Add the sentinel character.
	s.push_back(0);
}

/** Construct the suffix array. */
void FMIndex::buildSA(const vector<uint8_t>& s, vector<uint32_t>& sa)
{
	sa.resize(s.size());
	int status = saisxx(s.begin(), sa.begin(), (int)s.size(), 0x100);
	assert(status == 0);
	if (status != 0)
		abort();
}

/** Construct the Burrows–Wheeler transform. */
void FMIndex::buildBWT(const vector<uint8_t>& s,
		const vector<uint32_t>& sa, vector<uint64_t>& bwt)
{
	bwt.resize(s.size());
	for (size_t i = 0; i < s.size(); i++)
		bwt[i] = sa[i] > 0 ? s[sa[i] - 1] : s[s.size() - 1];
}

/** Sample the suffix array. */
void FMIndex::buildSampledSA(const vector<uint8_t>& s,
		const vector<uint32_t>& sa)
{
	for (size_t i = 0; i < s.size(); i++)
		if (i % m_sampleSA == 0)
			m_sampledSA.push_back(sa[i]);
}

/** Count the character frequency statistics. */
void FMIndex::calculateStatistics(const vector<uint8_t>& s)
{
	const unsigned UINT8_MAX = 255;
	vector<uint32_t> tmpCf(UINT8_MAX + 1);
	for (size_t i = 0; i < s.size(); i++)
		tmpCf[s[i]]++;
	m_cf.resize(UINT8_MAX + 1);
	m_cf[0] = 0;
	for (size_t i = 1; i <= UINT8_MAX; i++) {
		m_cf[i] = tmpCf[i-1];
		tmpCf[i] += tmpCf[i-1];
	}
}

/** Return the position of the specified suffix in the original
 * string.
 */
size_t FMIndex::locate(uint64_t i) const
{
	size_t bsize = m_wa.length();
	size_t j = i;
	size_t t = 0;
	while (j % m_sampleSA != 0) {
		unsigned c = m_wa.Lookup(j);
		j = m_cf[c] + m_wa.Rank(c, j + 1) - 1;
		t++;
	}
	if (m_sampledSA[j / m_sampleSA] + t >= bsize)
		return (size_t)m_sampledSA[j / m_sampleSA] + t - bsize;
	return (size_t)m_sampledSA[j / m_sampleSA] + t;
}

/** Store this index. */
void FMIndex::save(ostream& out)
{
	out.write((const char*)&m_sampleSA, sizeof m_sampleSA);
	out.write((const char*)&m_alphaSize, sizeof m_alphaSize);
	size_t cfSize = m_cf.size();
	out.write((const char*)&cfSize, sizeof cfSize);
	size_t sampledSASize = m_sampledSA.size();
	out.write((const char*)&sampledSASize, sizeof sampledSASize);
	out.write((const char*)&m_cf[0], sizeof m_cf[0]*m_cf.size());
	out.write((const char*)&m_sampledSA[0],
			m_sampledSA.size() * sizeof m_sampledSA[0]);

	vector<unsigned char> key;
	vector<uint8_t> val;
	for (unsigned i = 0; i < m_mapping.size(); ++i) {
		if (m_mapping[i] != UCHAR_MAX) {
			key.push_back(i);
			val.push_back(m_mapping[i]);
		}
	}
	uint32_t count = key.size();
	out.write((const char*)&count, sizeof count);
	out.write((const char*)&key[0], count * sizeof key[0]);
	out.write((const char*)&val[0], count * sizeof val[0]);

	m_wa.Save(out);
}

/** Load this index. */
void FMIndex::load(istream& in)
{
	in.read((char*)&m_sampleSA, sizeof m_sampleSA);
	in.read((char*)&m_alphaSize, sizeof m_alphaSize);
	size_t cfSize;
	in.read((char*)&cfSize, sizeof cfSize);
	size_t sampledSASize;
	in.read((char*)&sampledSASize, sizeof sampledSASize);
	m_cf.resize(cfSize);
	in.read((char*)&m_cf[0], cfSize * sizeof m_cf[0]);
	m_sampledSA.resize(sampledSASize);
	in.read((char*)&m_sampledSA[0],
			sampledSASize * sizeof m_sampledSA[0]);

	uint32_t count;
	in.read((char*)&count, sizeof count);

	vector<unsigned char> key;
	key.resize(count);
	in.read((char*)&key[0], count * sizeof key[0]);

	vector<uint8_t> val;
	val.resize(count);
	in.read((char*)&val[0], count * sizeof val[0]);

	m_mapping.clear();
	for (uint32_t i = 0; i < count; i++) {
		unsigned k = key[i];
		unsigned v = val[i];
		if (k >= m_mapping.size())
			m_mapping.resize(k + 1, UCHAR_MAX);
		assert(m_mapping[k] == UCHAR_MAX);
		m_mapping[k] = v;
	}

	m_wa.Load(in);
}

/** Build the character occurence table. */
void FMIndex::buildWaveletTree(const vector<uint64_t>& bwt)
{
	m_wa.Init(bwt);
}

/** Build the FM index. */
void FMIndex::buildFmIndex(const char* path, unsigned sampleSA)
{
	assert(sampleSA > 0);
	m_sampleSA = sampleSA;

	cerr << "start reading the input-file\n";
	vector<uint8_t> s;
	read(path, s);

	cerr << "alphabet size:" << (int)m_alphaSize << '\n';

	double sTime = clock();
	vector<uint32_t> sa;
	cerr << "build SA\n";
	buildSA(s, sa);

	cerr << "calculate statistics\n";
	calculateStatistics(s);

	vector<uint64_t> bwt;
	cerr << "build BWT\n";
	buildBWT(s, sa, bwt);

	cerr << "build WaveletTree\n";
	buildWaveletTree(bwt);

	cerr << "build sampledSA\n";
	buildSampledSA(s, sa);

	double eTime = clock();
	cerr << "cpu time: " << (eTime - sTime) / CLOCKS_PER_SEC << '\n';
}
