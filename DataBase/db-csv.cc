/**
 * Populate database and export to CSV
 */

#include "DB.h"

using namespace std;

static const char USAGE_MESSAGE[] =
"usage:\nabyss-db-csv inputDB.file table_name output.file\n";

static bool existFile(const char* filename)
{
	ifstream file(filename);
	return file;
}

template <typename D>
static bool existTable(
		D& db, const string& tablename)
{
	bool exists = false;
	dbVec V;
	V = db.readSqlToVec("select name from sqlite_master where type='table';");
	unsigned i = 0;
	while (i<V.size()) {
		if (tablename == V[i][0]) {
			exists = true;
			break;
		}
		i++;
	}
	return exists;
}

template  <typename D>
static string createCSV(
		D& db, const string& tablename, const char* csvname)
{
	ofstream csvfile(csvname);
	stringstream pStream, sStream, msg;
	string pstr, sstr;
	dbVec headers, values;

	pStream << "pragma table_info (" << tablename << ");";
	pstr = pStream.str();
	headers = db.readSqlToVec(pstr.c_str());
	for (unsigned i=0; i<headers.size(); i++) {
		csvfile << headers[i][1] << ((i == (headers.size()-1)) ? "" : ",");
	}
	csvfile << endl;

	sStream << "select * from " << tablename << ";";
	sstr = sStream.str();
	values = db.readSqlToVec(sstr.c_str());
	for (unsigned i=0; i<values.size(); i++) {
		for (unsigned j=0; j<values[i].size(); j++) {
			csvfile << values[i][j] << ((j == (values[i].size()-1)) ? "" : ",");
		}
		csvfile << endl;
	}
	msg << "Done! Check " << csvname << " now.\n";
	return msg.str();
}

int main(int argc, char** argv)
{
	if (argc != 4) {
		cout << USAGE_MESSAGE;
		exit(EXIT_FAILURE);
	} else if (!existFile (argv[1])) {
		cout << "Database file " << "'" << argv[1] << "' doesn't exist." << endl;
		cout << USAGE_MESSAGE;
		exit(EXIT_FAILURE);
	} else {
		DB db;
		init(db, argv[1]);
		string tablename(db.getProperTableName (argv[2]));
		if (!existTable (db, tablename)) {
			cout << "Table " << "'" << argv[2] << "' doesn't exist in " << argv[1] << "." << endl;
			cout << USAGE_MESSAGE;
			exit(EXIT_FAILURE);
		}
		cout << createCSV(db, tablename, argv[3]);
	}
}
