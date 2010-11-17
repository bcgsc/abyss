#ifndef IOUTIL_H
#define IOUTIL_H 1

#include <iostream>

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

#endif
