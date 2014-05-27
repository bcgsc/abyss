/**
 * A SQLite3 interface
 */
#ifndef DB_H
#define DB_H 1

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <vector>
#include <algorithm>
#include <sqlite3.h>

class DB {
 private:

	sqlite3* db;
	sqlite3_stmt* stmt;

	/** Provide some feedback when a SQL command is run. */
	static int callback (void*, int, char**, char**);

	void openDB (const char*, const int&);
	void closeDB ();

 public:
	typedef std::vector<std::vector<std::string> > db_mat;
	
	/** Verbose option inherited */
	int verbose_val;

	DB (const std::string& path, const int& v) {
		openDB (path.empty() ? "ABySS.db" : path.c_str(), v);
	}
	
	~DB () {
		closeDB ();
	}

	/** Manually toggle foreign_keys option on at each transaction.
	  * This limitation is specific to sqlite version ~3.
	  */
	std::string activeForeignKey (const std::string& s) {
		std::string s_pragma("pragma foreign_keys=on; ");
		return s_pragma += s;
	}
	
	bool query (const std::string& s) {
		char* errMsg = 0;	
		std::string new_s (activeForeignKey(s));	
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

	db_mat readSqlToVec (const std::string& s) {
		int cols, step;
		db_mat results;
		if (sqlite3_prepare (db, (const char*)s.c_str(), -1, &stmt, NULL) != SQLITE_OK)
			exit(EXIT_FAILURE);
		cols = sqlite3_column_count (stmt);
		while (true) { // for each row
			std::vector<std::string> temp;
			step = sqlite3_step (stmt);
			if (step == SQLITE_ROW) {
				for (int i=0; i<cols; i++)
					temp.push_back ((sqlite3_column_text (stmt, i) != NULL) ? std::string((const char*)sqlite3_column_text (stmt, i)) : "");
				results.push_back(temp);
			} else if (step == SQLITE_DONE) {
				break;
			} else {
				exit(EXIT_FAILURE);
			}
		}
		sqlite3_finalize (stmt);
		return results;
	}
	
	/** Return sql-friendly string that doesn't contain '-'. */
	std::string getProperTableName (const std::string& table) {
		std::string temp = table;
		replace(temp.begin(), temp.end(), '-', '_');
		return temp;
	}
};

#endif
