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
void FMIndex::read(const char* path, vector<T>& s)
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
	assert(s.size() < numeric_limits<size_type>::max());

	// Translate the alphabet.
	if (m_alphabet.empty())
		setAlphabet(s.begin(), s.end());
	transform(s.begin(), s.end(), s.begin(), Translate(*this));
	replace(s.begin(), s.end(), UCHAR_MAX, 0);
}

/** Construct the suffix array. */
void FMIndex::buildSA(const vector<T>& s, vector<size_type>& sa)
{
	sa.resize(s.size() + 1);
	sa[0] = s.size();
	int status = saisxx(s.begin(), sa.begin() + 1,
			(int)s.size(), 0x100);
	assert(status == 0);
	if (status != 0)
		abort();
}

/** Construct the Burrows–Wheeler transform. */
void FMIndex::buildBWT(const vector<T>& s,
		const vector<size_type>& sa, vector<T>& bwt)
{
	bwt.resize(sa.size());
	for (size_t i = 0; i < sa.size(); i++)
		bwt[i] = sa[i] == 0 ? SENTINEL() : s[sa[i] - 1];
}

/** Sample the suffix array. */
void FMIndex::buildSampledSA(const vector<size_type>& sa)
{
	for (size_t i = 0; i < sa.size(); i++)
		if (i % m_sampleSA == 0)
			m_sampledSA.push_back(sa[i]);
}

/** Count the character frequency statistics. */
void FMIndex::calculateStatistics(const vector<T>& s)
{
	const unsigned UINT8_MAX = 255;
	vector<size_type> tmpCf(UINT8_MAX + 1);
	for (size_t i = 0; i < s.size(); i++)
		tmpCf[s[i]]++;
	m_cf.resize(UINT8_MAX + 1);
	// The sentinel character occurs once.
	m_cf[0] = 1;
	for (size_t i = 0; i < UINT8_MAX; i++)
		m_cf[i + 1] = m_cf[i] + tmpCf[i];
}

/** Return the position of the specified suffix in the original
 * string.
 */
size_t FMIndex::locate(size_t i) const
{
	size_t bsize = m_occ.size();
	size_t j = i;
	size_t t = 0;
	while (j % m_sampleSA != 0) {
		T c = m_occ.at(j);
		j = c == SENTINEL() ? 0 : m_cf[c] + m_occ.rank(c, j + 1) - 1;
		t++;
	}
	if (m_sampledSA[j / m_sampleSA] + t >= bsize)
		return (size_t)m_sampledSA[j / m_sampleSA] + t - bsize;
	return (size_t)m_sampledSA[j / m_sampleSA] + t;
}

/** Build the FM index. */
void FMIndex::buildFmIndex(const char* path, unsigned sampleSA)
{
	assert(sampleSA > 0);
	m_sampleSA = sampleSA;

	cerr << "start reading the input-file\n";
	vector<T> s;
	read(path, s);

	cerr << "alphabet size:" << m_alphabet.size() << '\n';

	double sTime = clock();
	vector<size_type> sa;
	cerr << "build SA\n";
	buildSA(s, sa);

	cerr << "calculate statistics\n";
	calculateStatistics(s);

	vector<T> bwt;
	cerr << "build BWT\n";
	buildBWT(s, sa, bwt);

	cerr << "build the character occurence table\n";
	m_occ.assign(bwt);

	cerr << "build sampledSA\n";
	buildSampledSA(sa);

	double eTime = clock();
	cerr << "cpu time: " << (eTime - sTime) / CLOCKS_PER_SEC << '\n';
}
