#ifndef IOUTIL_H
#define IOUTIL_H 1

#include <cassert>
#include <cerrno>
#include <cstdlib>
#include <cstring> // for strerror
#include <fstream>
#include <iostream>
#include <limits> // for numeric_limits
#include <string>

/** Print an error message and exit if stream is not good. */
static inline void assert_good(const std::ios& stream,
		const std::string& path)
{
	if (!stream.good()) {
		std::cerr << "error: `" << path << "': "
			<< strerror(errno) << std::endl;
		exit(EXIT_FAILURE);
	}
}

/** Print an error message and exit if stream is not eof. */
static inline void assert_eof(std::istream& in,
		const std::string& path)
{
	if (!in.eof()) {
		in.clear();
		std::string s;
		std::getline(in, s);
		std::cerr << "error: `" << path << "': "
			"Expected end-of-file and saw `" << s << "'\n";
		exit(EXIT_FAILURE);
	}
}

/** This input stream manipulator skips the specified string. */
struct expect {
	const char* s;
	expect(const char* s) : s(s) { }
};

static inline std::istream& operator>>(std::istream& in, expect o)
{
	for (const char* p = o.s; *p != '\0'; ++p) {
		if (*p == ' ') {
			in >> std::ws;
		} else {
			char c = in.get();
			if (!in || c != *p) {
				std::cerr << "error: Expected `" << p
					<< "' and saw ";
				if (in) {
					std::cerr << '`' << c << "'\n";
					std::string s;
					if (getline(in, s) && !s.empty())
						std::cerr << "near: " << c << s << '\n';
				} else if (in.eof())
					std::cerr << "end-of-file\n";
				else
					std::cerr << "I/O error\n";
				exit(EXIT_FAILURE);
			}
		}
	}
	return in;
}

/** This input stream manipulator discards characters until reaching
 * the delimeter. */
struct Ignore {
	const char delim;
	size_t n;
	Ignore(const char delim,
			size_t n = std::numeric_limits<std::streamsize>::max())
		: delim(delim), n(n) { }
};

static inline std::istream& operator>>(std::istream& in, Ignore o)
{
	return in.ignore(o.n, o.delim);
}

/** Skip the specified character if it's next in the input stream. */
struct Skip {
	char c;
	Skip(const char c) : c(c) { }
};

static inline std::istream& operator>>(std::istream& in, Skip o)
{
	if (in.peek() == o.c)
		in.ignore(1);
	return in;
}

/** Read a file and store it in the specified vector. */
template <typename Vector>
static inline void readFile(const char* path, Vector& s)
{
	std::ifstream in(path);
	assert_good(in, path);
	in.seekg(0, std::ios::end);
	ssize_t n = in.tellg();
	assert(n > 0);
	s.resize(n);
	in.seekg(0, std::ios::beg);
	assert_good(in, path);
	char *p = reinterpret_cast<char*>(&s[0]);
#if __MACH__
	// Read 1 GB at a time. Reads of 2 GB or more fail.
	const ssize_t N = 1024 * 1024 * 1024;
	for (; n > N; n -= N, p += N) {
		in.read(p, N);
		assert_good(in, path);
		assert(in.gcount() == N);
	}
#endif
	in.read(p, n);
	assert_good(in, path);
	assert(in.gcount() == n);
}

/** Copy a file */
inline static void copyFile(const std::string& srcPath,
	const std::string& dstPath)
{
	assert(srcPath != dstPath);
	std::ifstream src(srcPath.c_str(), std::ios::binary);
	std::ofstream dst(dstPath.c_str(), std::ios::binary);
	dst << src.rdbuf();
	assert(dst);
}

#endif
