// DB.cc
#include "DB.h"
using namespace std;

int DB::callback (void* info, int argc, char** argv, char** colName)
{
	int i;
	cerr << (const char*)info << endl;
	for (i=0; i<argc; i++) {
		cerr << "%s = %s\n" << colName[i] << (argv[i] ? argv[i] : "NULL");
	}
	cerr << "\n";
	return 0;
}

void DB::openDB (const char* c, const int& v)
{
	verbose_val = v;
	if (sqlite3_open (c, &db)) {
		cerr << "Can't open DB.\n";
		exit(EXIT_FAILURE);
	}
	else {
		if (verbose_val > 0)
			cerr << "DB opened.\n";
	}
}

void DB::closeDB ()
{
	if (sqlite3_close (db)) {
		cerr << "Can't close DB.\n";
		exit(EXIT_FAILURE);
	} else {
		if (verbose_val > 0)
			cerr << "DB closed.\n";
	}
}
