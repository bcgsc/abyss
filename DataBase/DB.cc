#include "DB.h"

using namespace std;

void DB::closeDB()
{
	if (sqlite3_close(db)) {
		cerr << "[" << prog << "] Can't close DB.\n";
		exit(EXIT_FAILURE);
	} else {
		if (verbose_val >= 0 && exp != READ)
			cerr << "[" << prog << "] DB closed.\n";
	}
}

dbVec DB::readSqlToVec(const string& s)
{
	int cols, step;
	dbVec results;
	if (sqlite3_prepare(db, (const char*)s.c_str(), -1, &stmt, NULL) != SQLITE_OK)
		exit(EXIT_FAILURE);
	cols = sqlite3_column_count(stmt);
	while (true) {
		dbVars temp;
		step = sqlite3_step(stmt);
		if (step == SQLITE_ROW) {
			for (int i=0; i<cols; i++)
				temp.push_back(
					(sqlite3_column_text(stmt, i) != NULL) ? string((const char*)sqlite3_column_text(stmt, i)) : ""
				);
			results.push_back(temp);
		} else if (step == SQLITE_DONE) {
			break;
		} else {
			exit(EXIT_FAILURE);
		}
	}
	sqlite3_finalize(stmt);
	return results;
}

string DB::getProperTableName(const string& table)
{
	string temp = table;
	replace(temp.begin(), temp.end(), '-', '_');
	return temp;
}

void DB::createTables()
{
	stringstream sst;
	sst
		<< "\
create table if not exists Species (\
species_name text primary key);\
create table if not exists Strains (\
strain_name text primary key,species_name text,\
foreign key(species_name) references Species(species_name));\
create table if not exists Libraries (\
library_name text primary key,strain_name text,\
foreign key(strain_name) references Strains(strain_name));\
create table if not exists Run_pe (\
run_id text primary key,species_name text,strain_name text,library_name text,abyss_version text,\
pe_name text,pe_k integer,pe_lib text,\
time_start_run not null default (datetime(current_timestamp,'localtime')),stage integer default 0,\
foreign key(library_name) references Libraries(library_name));";

	if (query(sst.str()) && verbose_val > 2)
		cerr << sst.str() << endl;
}

void DB::insertToMetaTables(const dbVars& v)
{
	stringstream sst;
	sst
		<< "\
insert or ignore into species values('" << v[3] << "');\
insert or ignore into strains values('" << v[2] << "','" << v[3] << "');\
insert or ignore into libraries values('" << v[1] << "','" << v[2] << "');\
insert into run_pe(run_id,species_name,strain_name,library_name,abyss_version,pe_name,pe_k,pe_lib) \
values('" << v[0] << "','" << v[3] << "','" << v[2] << "','" << v[1] << "','" << VERSION << "','" << v[4] <<
 "','" << v[5] << "','" << v[6] << "');";

	if (query(sst.str()) && verbose_val > 1)
		cerr << sst.str() << endl;
}

string DB::initializeRun()
{
	string id("");
	createTables();
	ifstream ifile("db.txt");
	if (ifile) {
		dbVars iV;
		string eachLine;
		while (getline(ifile, eachLine)) iV.push_back(eachLine);
		if (iV.size() == 7) {
			insertToMetaTables(iV);
			ofstream ofile("db.txt");
			ofile << iV[0] << "\n";
		}
		stringstream uStream;
		uStream
			<< "update Run_pe set stage=stage+1 where run_id = '" << iV[0] << "';";
		if (query (uStream.str()) && verbose_val > 2)
			cerr << uStream.str() << endl;
		id = iV[0];
	}
	return id;
}

string DB::getPath(const string& program)
{
	stringstream echoStream, pathStream;
	echoStream
		<< "echo `which " << program << "` > temp.txt";
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

bool DB::definePeVars(const string& id)
{
	bool defined = false;
	stringstream select_cmd;
	select_cmd
		<< "select species_name, strain_name, library_name from Run_pe where run_id = '" << id << "';";
	dbVec v;
	v = readSqlToVec(select_cmd.str());
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

void DB::assemblyStatsToDb()
{
	string temp = getProperTableName(prog);
	stringstream sqlStream, mapKeys, mapValues;
	string id = initializeRun();

	if (temp == "ABYSS" || temp == "ABYSS_P")
		temp = "ABYSS_assembly";

	sqlStream
		<< "\
create table if not exists " << temp << " (\
run_id text default null, run_stage integer default null, species_name text, strain_name text, library_name text, \
exec_path text, command_line text, time_finish_" << temp << " not null default (datetime(current_timestamp,'localtime')), ";

	mapKeys
		<< "run_id, run_stage, species_name, strain_name, library_name, exec_path, command_line, ";

	if (id.length() == 0)
		mapValues
			<< "null,null,";
	else
		mapValues
			<< "'" << id << "',(select stage from Run_pe where run_id = '" << id << "'),";

	if ((id.length() > 0) && definePeVars(id))
		mapValues
			<< "'" << peVars[0] << "', '" << peVars[1] << "', '" << peVars[2] << "', '";
	else
		mapValues
			<< "'" << initVars[2] << "', '" << initVars[1] << "', '" << initVars[0] << "', '";

	mapValues
		<< getPath(prog)
		<< "', '" << cmd << "', ";

	while (!statMap.empty()) {
		mapKeys
			<< statMap.getFirst(statMap.begin())
			<< (statMap.size() == 1 ? "" : ", ");
		mapValues
			<< statMap.getSecond(statMap.begin())
			<< (statMap.size() == 1 ? "" : ", ");
		sqlStream
			<< statMap.getFirst(statMap.begin())
			<< (statMap.size() == 1 ? " int, foreign key(run_id) references Run_pe(run_id));" : " int, ");
		statMap.erase(statMap.begin());
	}

	if (query(sqlStream.str()) && verbose_val > 1)
		cerr << sqlStream.str() << endl;

	sqlStream.str("");
	sqlStream
		<< "insert into "<< temp << " (" << mapKeys.str() << ") values (" << mapValues.str() << ");";

	if (query(sqlStream.str()) && verbose_val > 0)
		cerr << sqlStream.str() << endl;
}
