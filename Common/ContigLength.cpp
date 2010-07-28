#include "ContigLength.h"
#include "ContigID.h"
#include <cassert>
#include <cerrno>
#include <cstdlib>
#include <cstring> // for strerror
#include <fstream>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

using namespace std;

namespace opt {
	extern unsigned k;
}

static void assert_open(ifstream& f, const string& p)
{
	if (f.is_open())
		return;
	cerr << p << ": " << strerror(errno) << endl;
	exit(EXIT_FAILURE);
}

/** Read contig lengths. */
vector<unsigned> readContigLengths(istream& in)
{
	assert(in);
	assert(ContigID::empty());
	vector<unsigned> lengths;
	ContigID id;
	unsigned len;
	while (in >> id >> len) {
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
	assert_open(in, path);
	return readContigLengths(in);
}
