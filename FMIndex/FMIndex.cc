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
