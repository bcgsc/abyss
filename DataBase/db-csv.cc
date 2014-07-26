/**
 * Populate SQLite database table(s) and export to CSV
 */

#include "DB.h"
#include <string.h>

using namespace std;

static const char USAGE_MESSAGE[] =
"usage:\
\nabyss-db-csv SQLite_repository table_name(s)\
\nabyss-db-csv SQLite_repository --all\n";

static const char TABLE_LIST[] =
"select name from sqlite_master where type='table';";

typedef vector<string> vs;

static bool existFile(const char* f)
{
	ifstream file(f);
	return file;
}

template <typename D>
static bool existTable(
		D& db, const string& t)
{
	dbVec v = db.readSqlToVec(TABLE_LIST);
	unsigned i = 0;

	while (i<v.size()) {
		if (t == v[i][0])
			return true;
		i++;
	}
	return false;
}

template  <typename D>
static string createCSV(
		D& db, const string& t, const string& csv)
{
	ofstream csvf(csv.c_str());
	stringstream ps, ss, msg;
	string pstr, sstr;
	dbVec head, val;

	ps << "pragma table_info (" << t << ");";
	pstr = ps.str();
	head = db.readSqlToVec(pstr.c_str());

	for (unsigned i=0; i<head.size(); i++)
		csvf << head[i][1] << ((i == (head.size()-1)) ? "" : ",");
	csvf << endl;

	ss << "select * from " << t << ";";
	sstr = ss.str();
	val = db.readSqlToVec(sstr.c_str());

	for (unsigned i=0; i<val.size(); i++) {
		for (unsigned j=0; j<val[i].size(); j++)
			csvf << val[i][j] << ((j == (val[i].size()-1)) ? "" : ",");
		csvf << endl;
	}
	msg << csv << " has been created.\n";
	return msg.str();
}

template <typename D>
static void populateOneTable(
		D& db, const string& repo, const string& t)
{
	string tt(db.getProperTableName(t));
	stringstream csv;
	csv << "db." << tt << ".csv";

	if (!existTable(db, tt)) {
		cout << "Table " << "'" <<
			tt << "' doesn't exist in " << repo << "." << endl;
		cout << USAGE_MESSAGE;
		exit(EXIT_FAILURE);
	}
	cout << createCSV(db, tt, csv.str());
}

template <typename D>
static void populateAll(
		D& db, const string& repo)
{
	dbVec v = db.readSqlToVec(TABLE_LIST);

	for (unsigned i=0; i<v.size(); i++)
		populateOneTable(db, repo, v[i][0]);
}

int main(int argc, char** argv)
{
	if (argc < 3) {
		cout << USAGE_MESSAGE;
		exit(EXIT_FAILURE);
	} else if (!existFile(argv[1])) {
		cout << "Database file " << "'" <<
			argv[1] << "' doesn't exist." << endl;
		cout << USAGE_MESSAGE;
		exit(EXIT_FAILURE);
	} else {
		DB db;
		init(db, argv[1]);

		if (argc == 3 && strcmp(argv[2], "--all") == 0) {
			populateAll(db, argv[1]);
			exit(EXIT_SUCCESS);
		}

		vs t(argc-2);
		char** last = argv + argc;
		copy(argv+2, last, t.begin());

		for (unsigned i=0; i<t.size(); i++)
			populateOneTable(db, argv[1], t[i]);
	}
}
