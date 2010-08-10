#include "FastaReader.h"
#include "DataLayer/Options.h"
#include <algorithm>
#include <cassert>
#include <cctype>
#include <cerrno>
#include <cstdlib>
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

	/** minimum quality threshold */
	int qualityThreshold;

	/** quality offset, usually 33 or 64 */
	int qualityOffset;
}

/** Output an error message. */
ostream& FastaReader::die()
{
	return cerr << m_path << ':' << m_line << ": error: ";
}

static void assert_open(ifstream& f, const string& p)
{
	if (f.is_open())
		return;
	cerr << p << ": error: " << strerror(errno) << endl;
	exit(EXIT_FAILURE);
}

FastaReader::FastaReader(const char* path, int flags)
	: m_path(path), m_fin(path),
	m_in(strcmp(path, "-") == 0 ? cin : m_fin),
	m_flags(flags), m_line(0), m_unchaste(0)
{
	if (strcmp(path, "-") != 0)
		assert_open(m_fin, path);
	if (m_in.peek() == EOF)
		cerr << m_path << ':' << m_line << ": warning: "
			"file is empty\n";
}

/** Return whether this read passed the chastity filter. */
bool FastaReader::isChaste(const string& s, const string& line)
{
	if (s == "1" || s == "Y") {
		return true;
	} else if (s == "0" || s == "N") {
		return false;
	} else {
		die() << "chastity filter should be one of 0, 1, N or Y\n"
			"and saw `" << s << "' near\n" << line << endl;
		exit(EXIT_FAILURE);
	}
}

/** Check that the seqeuence and quality agree in length. */
void FastaReader::checkSeqQual(const string& s, const string& q)
{
	if (s.length() != q.length()) {
		die() << "sequence and quality must be the same length near\n"
			<< s << '\n' << q << endl;
		exit(EXIT_FAILURE);
	}
}

/** Read a single record. */
Sequence FastaReader::read(string& id, string& comment,
		char& anchor, string& q)
{
next_record:
	// Discard comments.
	while (m_in.peek() == '#')
		m_in.ignore(numeric_limits<streamsize>::max(), '\n');

	signed char recordType = m_in.peek();
	Sequence s;

	unsigned qualityOffset = 0;
	if (recordType == EOF) {
		m_in.get();
		return s;
	} else if (recordType == '>' || recordType == '@') {
		// Read the header.
		string header;
		getline(header);
		istringstream headerStream(header);
		headerStream >> recordType >> id >> ws;
		std::getline(headerStream, comment);

		// Ignore SAM headers.
		if (id.length() == 2 && isupper(id[0]) && isupper(id[1])
				&& comment.length() > 2 && comment[2] == ':')
			goto next_record;

		getline(s);

		if (recordType == '@') {
			char c = m_in.get();
			if (c != '+') {
				string line;
				getline(line);
				die() << "expected `+' and saw `" << c << "' near\n"
					<< c << line << "\n^\n";
				exit(EXIT_FAILURE);
			}
			m_in.ignore(numeric_limits<streamsize>::max(), '\n');
			getline(q);
			checkSeqQual(s, q);
		} else
			q.clear();

		if (opt::trimMasked) {
			// Removed masked (lower case) sequence at the beginning
			// and end of the read.
			size_t trimFront = s.find_first_not_of("acgtn");
			size_t trimBack = s.find_last_not_of("acgtn") + 1;
			s.erase(trimBack);
			s.erase(0, trimFront);
			if (!q.empty()) {
				q.erase(trimBack);
				q.erase(0, trimFront);
			}
		}
		if (flagFoldCase())
			transform(s.begin(), s.end(), s.begin(), ::toupper);

		if (s.length() > 2 && isalpha(s[0]) && isdigit(s[1])) {
			// The first character is the primer base. The second
			// character is the dibase read of the primer and the
			// first base of the sample, which is not part of the
			// assembly.
			anchor = colourToNucleotideSpace(s[0], s[1]);
			s.erase(0, 2);
		}

		qualityOffset = 33;
	} else {
		string line;
		vector<string> fields;
		fields.reserve(22);
		getline(line);
		istringstream in(line);
		string field;
		while (std::getline(in, field, '\t'))
			fields.push_back(field);

		if (fields.size() >= 11
				&& fields[9].length() == fields[10].length()) {
			// SAM
			unsigned flags = strtoul(fields[1].c_str(), NULL, 0);
			if (flags & 0x100) // FSECONDARY
				goto next_record;
			if (opt::chastityFilter && (flags & 0x200)) { // FQCFAIL
				m_unchaste++;
				goto next_record;
			}
			id = fields[0];
			switch (flags & 0xc1) { // FPAIRED|FREAD1|FREAD2
			  case 0: case 1: break; // FPAIRED
			  case 0x41: id += "/1"; break; // FPAIRED|FREAD1
			  case 0x81: id += "/2"; break; // FPAIRED|FREAD2
			  default:
				die() << "invalid flags: `" << id << "' near"
					<< line << endl;
				exit(EXIT_FAILURE);
			}
			s = fields[9];
			q = fields[10];
			if (s == "*")
				s.clear();
			if (q == "*")
				q.clear();
			if (flags & 0x10) { // FREVERSE
				s = reverseComplement(s);
				reverse(q.begin(), q.end());
			}
			comment = fields[1];
			qualityOffset = 33;
			if (!q.empty())
				checkSeqQual(s, q);
		} else if (fields.size() == 11 || fields.size() == 22) {
			// qseq or export
			if (opt::chastityFilter
					&& !isChaste(fields.back(), line)) {
				m_unchaste++;
				goto next_record;
			}

			ostringstream o;
			o << fields[0];
			for (int i = 1; i < 6; i++)
				if (!fields[i].empty())
					o << ':' << fields[i];
			if (!fields[6].empty() && fields[6] != "0")
				o << '#' << fields[6];
			// The reverse read is typically the second read, but is
			// the third read of an indexed run.
			o << '/' << (fields[7] == "3" ? "2" : fields[7]);
			id = o.str();
			s = fields[8];
			q = fields[9];
			comment = fields.back();
			qualityOffset = 64;
			checkSeqQual(s, q);
		} else {
			die() << "Expected either `>' or `@' or 11 fields\n"
					"and saw `" << recordType << "' and "
					<< fields.size() << " fields near\n"
					<< line << endl;
			exit(EXIT_FAILURE);
		}
	}

	if (opt::qualityOffset > 0)
		qualityOffset = opt::qualityOffset;

	if (opt::qualityThreshold > 0 && !q.empty()) {
		assert(s.length() == q.length());
		static const char ASCII[] =" !\"#$%&'()*+,-./0123456789"
			":;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefgh";
		assert(qualityOffset > (unsigned)ASCII[0]);
		const char* goodQual = ASCII + (qualityOffset - ASCII[0])
			+ opt::qualityThreshold;

		size_t trimFront = q.find_first_of(goodQual);
		size_t trimBack = q.find_last_of(goodQual) + 1;
		if (trimFront > 0 || trimBack < q.length()) {
			s.erase(trimBack);
			s.erase(0, trimFront);
		}
	}

	assert(qualityOffset >= 33);
	if (flagConvertQual() && qualityOffset != 33) {
		// Convert to standard quality (ASCII 33).
		for (string::iterator it = q.begin(); it != q.end(); ++it) {
			int x = *it - qualityOffset;
			if (x < -5 || x > 40) {
				die() << "quality " << x
					<< " is out of range -5 <= q <= 40 near\n"
					<< q << '\n'
					<< string(it - q.begin(), ' ') << "^\n";
				exit(EXIT_FAILURE);
			}
			*it = 33 + max(0, x);
		}
	}

	return s;
}
