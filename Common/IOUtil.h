#ifndef IOUTIL_H
#define IOUTIL_H 1

#include <iostream>
#include <limits> // for numeric_limits

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
				std::cerr << "error: Expected `" << *p
					<< "' and saw ";
				if (in)
					std::cerr << '`' << c << "'\n";
				else if (in.eof())
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
struct ignore {
	const char delim;
	size_t n;
	ignore(const char delim,
			size_t n = std::numeric_limits<std::streamsize>::max())
		: delim(delim), n(n) { }
};

static inline std::istream& operator>>(std::istream& in, ignore o)
{
	return in.ignore(o.n, o.delim);
}

#endif
