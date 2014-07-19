/**
 * A SQLite3 interface
 */

#ifndef DB_H
#define DB_H 1

#include "Common/InsOrderedMap.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <config.h>
#include <stdlib.h>
#include <sqlite3.h>

typedef std::vector<std::vector<std::string> > dbVec; // for output
typedef InsOrderedMap<std::string,int> dbMap; // for input
typedef std::vector<std::string> dbVars;

class DB {
 private:
	dbMap statMap;
	dbVars initVars, peVars;
	sqlite3* db;
	sqlite3_stmt* stmt;
	std::string prog, cmd;
	int exp;

	// Verbosity inherited from the equivalent abyss-pe option.
	int verbose_val;

	void openDB (const char* c, const int& v) {
		verbose_val = v;
		if (sqlite3_open (c, &db)) {
			std::cerr << "Can't open DB.\n";
			exit(EXIT_FAILURE);
		} else {
			if (verbose_val > 0)
				std::cerr << "DB opened.\n";
		}
	}

	static int callback (void* info, int argc, char** argv, char** colName) {
		int i;
		std::cerr << (const char*)info << std::endl;
		for (i=0; i<argc; i++) {
			std::cerr << "%s = %s\n" << colName[i] << (argv[i] ? argv[i] : "NULL");
		}
		std::cerr << "\n";
		return 0;
	}

	void closeDB ();
	void createTables ();
	void insertToMetaTables (const dbVars&);
	bool isRun ();
	std::string getPath (const std::string&);
	bool definePeVars ();
	void assemblyStatsToDb ();

 public:
	DB () {
		initVars.assign (3,"");
		peVars.assign (3,"");
	}

	~DB () {
		if (exp == 0)
			assemblyStatsToDb ();
		closeDB ();
	}

	friend void init (DB& d, const std::string& path) {
		d.openDB (path.c_str(), 0);
		d.exp = 1;
	}

	friend void init (DB& d, const std::string& path, const int& v, const std::string& program, const std::string& command, const dbVars& vars) {
		// If destination is not specified, create 'ABySS.db' by default.
		d.openDB (path.empty() ? "ABySS.db" : path.c_str(), v);
		d.prog = program;
		d.cmd = command;
		d.initVars = vars;
		d.exp = 0;
	}

	std::string activateForeignKey (const std::string& s) {
		std::string s_pragma("pragma foreign_keys=on; ");
		return s_pragma += s;
	}

	bool query (const std::string& s) {
		char* errMsg = 0;
		std::string new_s (activateForeignKey(s));
		const char* statement = new_s.c_str();
		int rc = sqlite3_exec (db, statement, callback, 0, &errMsg);
		if (rc != SQLITE_OK) {
			std::cerr << "SQL error: " << errMsg << std::endl;
			sqlite3_free (errMsg);
			exit(EXIT_FAILURE);
		} else {
			return true;
		}
	}

	friend void addToDb (DB& d, const std::string& key, const int& value) {
		d.statMap.push_back (key, value);
	}

	friend void addToDb (DB& d, dbMap m) {
		d.statMap.insert (m.getAC());
	}

	dbVec readSqlToVec (const std::string&);
	std::string getProperTableName (const std::string&); 
};

#endif
