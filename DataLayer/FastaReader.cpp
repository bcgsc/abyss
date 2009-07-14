#include "FastaReader.h"
#include "Options.h"
#include <algorithm>
#include <cassert>
#include <cctype>
#include <cerrno>
#include <cstring>
#include <sstream>
#include <iostream>

using namespace std;

static void assert_open(ifstream& f, const string& p)
{
	if (f.is_open())
		return;
	cerr << p << ": " << strerror(errno) << endl;
	exit(EXIT_FAILURE);
}

FastaReader::FastaReader(const char* path)
	: m_inPath(path), m_inFile(path),
	m_fileHandle(strcmp(path, "-") == 0 ? cin : m_inFile),
	m_nonacgt(0)
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

Sequence FastaReader::ReadSequence(string& id)
{
	// Discard comments.
	while (m_fileHandle.peek() == '#') {
		m_fileHandle.ignore(numeric_limits<streamsize>::max(), '\n');
		if (m_fileHandle.peek() == EOF) {
			fputs("error: file ends in comments\n", stderr);
			assert(false);
		}
	}

	char recordType = m_fileHandle.peek();
	Sequence s;

	if (recordType == '>' || recordType == '@') {
		// Read the header.
		m_fileHandle >> recordType >> id;
		m_fileHandle.ignore(numeric_limits<streamsize>::max(), '\n');

		getline(m_fileHandle, s);
		transform(s.begin(), s.end(), s.begin(), ::toupper);

		assert(s.length() > 2);
		if (isalpha(s[0]) && isdigit(s[1])) {
			// The first character is the primer base. The second
			// character is the dibase read of the primer and the first
			// base of the sample, which is not part of the assembly.
			s = s.substr(2);
		}

		if (recordType == '@') {
			// Discard the quality values.
			char c;
			m_fileHandle >> c;
			assert(c == '+');
			m_fileHandle.ignore(numeric_limits<streamsize>::max(), '\n');
			m_fileHandle.ignore(numeric_limits<streamsize>::max(), '\n');
		}
	} else {
		string line;
		vector<string> fields;
		fields.reserve(11);
		getline(m_fileHandle, line);
		istringstream in(line);
		string field;
		while (in >> field)
			fields.push_back(field);

		if (fields.size() == 11) {
			id = fields[0];
			for (int i = 1; i < 6; i++) {
				id.append("_");
				id.append(fields[i]);
			}
			id.append("/");
			id.append(fields[7]);
			s = fields[8];
		} else {
			fprintf(stderr, "error: `%s' is an unknown format\n"
					"Expected either `>' or `@' or 11 fields\n"
					"and saw `%c' and %d fields\n",
					m_inPath, recordType, fields.size());
			exit(EXIT_FAILURE);
		}
	}

	return s;
}

// Read in a group of sequences and return whether there are sequences remaining
bool FastaReader::ReadSequences(SequenceVector& outseqs)
{
	if (!isGood())
		return false;
	Sequence seq = ReadSequence();
	size_t pos = seq.find_first_not_of("ACGT0123");
	if (pos == std::string::npos) {
		outseqs.push_back(seq);
	} else {
		if (opt::verbose > 4)
			fprintf(stderr,
					"warning: discarded sequence containing `%c'\n",
					seq[pos]);
		m_nonacgt++;
	}
	return true;
}
