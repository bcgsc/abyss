#include "DB.h"
using namespace std;

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

dbVec DB::readSqlToVec (const string& s)
{
	int cols, step;
	dbVec results;
	if (sqlite3_prepare (db, (const char*)s.c_str(), -1, &stmt, NULL) != SQLITE_OK)
		exit(EXIT_FAILURE);
	cols = sqlite3_column_count (stmt);
	while (true) {
		dbVars temp;
		step = sqlite3_step (stmt);
		if (step == SQLITE_ROW) {
			for (int i=0; i<cols; i++)
				temp.push_back ((sqlite3_column_text (stmt, i) != NULL) ? string((const char*)sqlite3_column_text (stmt, i)) : "");
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

string DB::getProperTableName (const string& table)
{
	string temp = table;
	replace(temp.begin(), temp.end(), '-', '_');
	return temp;
}

void DB::createTables ()
{
	stringstream sst;
	sst << "\
create table if not exists Species (species_name text primary key);\
create table if not exists Strains (strain_name text primary key,species_name text,\
foreign key(species_name) references Species(species_name));\
create table if not exists Libraries (library_name text primary key,strain_name text,\
foreign key(strain_name) references Strains(strain_name));\
create table if not exists Run_pe (run_id integer primary key,species_name text,strain_name,library_name text,abyss_version text,time_start_run not null default (datetime(current_timestamp,'localtime')),stage integer default 0,\
foreign key(library_name) references Libraries(library_name));";

	if (query (sst.str()) && verbose_val > 1)
		cerr << sst.str() << endl;
}

void DB::insertToMetaTables (const dbVars& v)
{
	stringstream sst;
	sst << "\
insert or ignore into species values('" << v[2] << "');\
insert or ignore into strains values('" << v[1] << "','" << v[2] << "');\
insert or ignore into libraries values('" << v[0] << "','" << v[1] << "');\
insert into run_pe(run_id,species_name,strain_name,library_name,abyss_version) values(null,'" << v[2] << "','" << v[1] << "','" << v[0] << "','" << VERSION << "');";

	if (query (sst.str()) && verbose_val > 1)
		cerr << sst.str() << endl;
}

bool DB::isRun ()
{
	bool isrun = false;
	createTables ();
	ifstream ifile("db.txt");
	if (ifile) {
		dbVars iV;
		string eachLine;
		while (getline (ifile, eachLine)) iV.push_back(eachLine);
		if (iV.size() == 3) {
			insertToMetaTables (iV);
			ofstream ofile("db.txt");
			ofile << "done\n";
		}
		stringstream uStream;
		uStream << "update Run_pe set stage=stage+1 where run_id = (select max(run_id) from Run_pe);";
		if (query (uStream.str()) && verbose_val >2)
			cerr << uStream.str() << endl;
		isrun = true;
	}
	return isrun;
}

string DB::getPath (const string& program)
{
	stringstream echoStream, pathStream;
	echoStream << "echo `which " << program << "` > temp.txt";
	string echoStr = echoStream.str();
	system(echoStr.c_str());
	ifstream ifile("temp.txt");
	if (ifile) {
		string line;
		while (getline (ifile, line)) pathStream << line;
	}
	system("rm -f temp.txt");
	return pathStream.str();
}

bool DB::definePeVars ()
{
	bool defined = false;
	dbVec v;
	v = readSqlToVec ("select species_name, strain_name, library_name from Run_pe where run_id = (select max(run_id) from Run_pe);");
	if (v[0].size() == peVars.size()) {
		unsigned i = 0;
		while (i<v[0].size()) {
			peVars[i] = v[0][i];
			i++;
		}
		defined = true;
	}
	return defined;
}

void DB::assemblyStatsToDb ()
{
	string temp = getProperTableName(prog);
	stringstream sqlStream, mapKeys, mapValues;
	bool isrun = isRun();

	sqlStream << "\
create table if not exists " << temp << " (\
run_id integer default null, run_stage integer default null, species_name text, strain_name text, library_name text, \
exec_path text, command_line text, time_finish_" << temp << " not null default (datetime(current_timestamp,'localtime')), ";

	mapKeys << "run_id, run_stage, species_name, strain_name, library_name, exec_path, command_line, ";
	mapValues << (isrun ? "(select max(run_id) from Run_pe), (select stage from Run_pe where run_id = (select max(run_id) from Run_pe)), " : "null, null, ");
	if (isrun && definePeVars())
		mapValues  << "'" << peVars[2] << "', '" << peVars[1] << "', '" << peVars[0] << "', '";
	else
		mapValues  << "'" << initVars[2] << "', '" << initVars[1] << "', '" << initVars[0] << "', '";
	mapValues << getPath(prog) << "', '" << cmd << "', ";

	while (!statMap.empty()) {
		mapKeys << statMap.begin()->first << (statMap.size() == 1 ? "" : ", ");
		mapValues << statMap.begin()->second << (statMap.size() == 1 ? "" : ", ");
		sqlStream << statMap.begin()->first << (statMap.size() == 1 ? " int, foreign key(run_id) references Run_pe(run_id));" : " int, ");
		statMap.erase(statMap.begin());
	}

	if (query (sqlStream.str()) && verbose_val > 1)
		cerr << sqlStream.str() << endl;

	sqlStream.str("");
	sqlStream << "insert into " << temp << " (" << mapKeys.str() << ") values (" << mapValues.str() << ");";

	if (query (sqlStream.str()) && verbose_val > 1)
		cerr << sqlStream.str() << endl;
}
