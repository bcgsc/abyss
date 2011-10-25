#include "ContigLength.h"
#include "ContigID.h"
#include "IOUtil.h"
#include <cassert>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

using namespace std;

namespace opt {
	extern unsigned k;
}

/** Read contig lengths. */
vector<unsigned> readContigLengths(istream& in)
{
	assert(in);
	assert(ContigID::empty());
	vector<unsigned> lengths;
	string s;
	unsigned len;
	while (in >> s >> len) {
		ContigID id = ContigID::insert(s);
		in.ignore(numeric_limits<streamsize>::max(), '\n');
		assert(len >= opt::k);
		assert(id == lengths.size());
		lengths.push_back(len - opt::k + 1);
	}
	assert(in.eof());
	assert(!lengths.empty());
	ContigID::lock();
	return lengths;
}

/** Read contig lengths. */
vector<unsigned> readContigLengths(const string& path)
{
	ifstream in(path.c_str());
	assert_good(in, path);
	return readContigLengths(in);
}
