#include "FastaReader.h"
#include "DataLayer/Options.h"
#include "Log.h"
#include <algorithm>
#include <cassert>
#include <cctype>
#include <cerrno>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <limits> // for numeric_limits
#include <sstream>
#include <vector>

using namespace std;

namespace opt {
	/** Discard reads that failed the chastity filter. */
	int chastityFilter = 1;

	/** Trim masked (lower case) characters from the ends of
	 * sequences.
	 */
	int trimMasked = 1;
}

static void assert_open(ifstream& f, const string& p)
{
	if (f.is_open())
		return;
	cerr << p << ": " << strerror(errno) << endl;
	exit(EXIT_FAILURE);
}

FastaReader::FastaReader(const char* path, bool discardN)
	: m_inPath(path), m_inFile(path),
	m_fileHandle(strcmp(path, "-") == 0 ? cin : m_inFile),
	m_discardN(discardN),
	m_unchaste(0), m_nonacgt(0)
{
	if (strcmp(path, "-") != 0)
		assert_open(m_inFile, path);
	if (m_fileHandle.peek() == EOF)
		fprintf(stderr, "warning: `%s' is empty\n", path);
}

FastaReader::~FastaReader()
{
	m_inFile.close();
}

static bool isChaste(const string &s)
{
	if (s == "1" || s == "Y") {
		return true;
	} else if (s == "0" || s == "N") {
		return false;
	} else {
		cerr << "error: chastity filter should be either "
			<< "0, 1, N or Y and saw `" << s << "'\n";
		exit(EXIT_FAILURE);
	}
}

Sequence FastaReader::read(string& id, string& comment, char& anchor)
{
next_record:
	// Discard comments.
	while (m_fileHandle.peek() == '#')
		m_fileHandle.ignore(numeric_limits<streamsize>::max(), '\n');

	signed char recordType = m_fileHandle.peek();
	Sequence s;

	if (recordType == EOF) {
		m_fileHandle.get();
		return s;
	} else if (recordType == '>' || recordType == '@') {
		// Read the header.
		string header;
		getline(m_fileHandle, header);
		stringstream headerStream(header);
		headerStream >> recordType >> id >> ws;
		getline(headerStream, comment);

		getline(m_fileHandle, s);
		assert(s.length() > 0);

		if (opt::trimMasked) {
			// Removed masked (lower case) sequence at the beginning
			// and end of the read.
			s.erase(s.find_last_not_of("acgtn") + 1);
			s.erase(0, s.find_first_not_of("acgtn"));
		}
		transform(s.begin(), s.end(), s.begin(), ::toupper);

		if (s.length() > 2 && isalpha(s[0]) && isdigit(s[1])) {
			// The first character is the primer base. The second
			// character is the dibase read of the primer and the
			// first base of the sample, which is not part of the
			// assembly.
			anchor = colourToNucleotideSpace(s[0], s[1]);
			s.erase(0, 2);
		}

		if (recordType == '@') {
			// Discard the quality values.
			char c = m_fileHandle.get();
			assert(c == '+');
			(void)c;
			m_fileHandle.ignore(numeric_limits<streamsize>::max(), '\n');
			m_fileHandle.ignore(numeric_limits<streamsize>::max(), '\n');
		}
	} else {
		string line;
		vector<string> fields;
		fields.reserve(22);
		getline(m_fileHandle, line);
		istringstream in(line);
		string field;
		while (getline(in, field, '\t'))
			fields.push_back(field);

		if (fields.size() == 11 || fields.size() == 22) {
			if (opt::chastityFilter && !isChaste(fields.back())) {
				m_unchaste++;
				goto next_record;
			}

			ostringstream o;
			o << fields[0];
			for (int i = 1; i < 6; i++)
				o << '_' << fields[i];
			o << '/' << fields[7];
			id = o.str();
			s = fields[8];
		} else {
			fprintf(stderr, "error: `%s' is an unknown format\n"
					"Expected either `>' or `@' or 11 fields\n"
					"and saw `%c' and %zu fields\n",
					m_inPath, recordType, fields.size());
			exit(EXIT_FAILURE);
		}
	}

	if (m_discardN) {
		size_t pos = s.find_first_not_of("ACGT0123");
		if (pos != string::npos) {
			logger(5) << "warning: discarded sequence containing `"
				<< s[pos] << "'\n";
			m_nonacgt++;
			goto next_record;
		}
	}

	return s;
}
