#ifndef STATQUERIES_H
#define STATQUERIES_H 1

#include "DB.h"
#include <sstream>
#include <config.h>
#include <map>
#include <stdlib.h>

static std::map<std::string,float> dbMap;

/** Create sql tables of species, libraries and abyss versions. */
template <typename D>
static inline void createTables (D& db)
{
	std::stringstream createStream;
	createStream << "\
create table if not exists Species (species_name text primary key);\
create table if not exists Strains (strain_name text primary key,species_name text,\
foreign key(species_name) references Species(species_name));\
create table if not exists Libraries (library_name text primary key,strain_name text,\
foreign key(strain_name) references Strains(strain_name));\
create table if not exists Abyss_pe (abyss_version text primary key);\
create table if not exists Run_pe (run_id integer primary key,library_name text,abyss_version text,time_start_run not null default (datetime(current_timestamp,'localtime')),stage integer default 0,\
foreign key(library_name) references Libraries(library_name),\
foreign key(abyss_version) references Abyss_pe(abyss_version));";

	if (db.query (createStream.str()) && db.verbose_val > 1)
		std::cerr << createStream.str() << std::endl;
}

/** Tell if current operation is run by abyss-pe. */
template <typename D>
static inline bool isRun (D& db)
{
	bool isrun = false;
	createTables (db);
	std::ifstream ifile("db.txt");
	if (ifile) {
		std::vector<std::string> iV;
		std::string eachLine;
		while (getline (ifile, eachLine)) iV.push_back(eachLine);

		if (iV.size() == 3) { // when species, strain, library are set
			std::stringstream iStream;
			iStream << "\
insert or ignore into species values('" << iV[0] << "');\
insert or ignore into strains values('" << iV[1] << "','" << iV[0] << "');\
insert or ignore into libraries values('" << iV[2] << "','" << iV[1] << "');\
insert or ignore into abyss_pe values('" << VERSION << "');\
insert into run_pe(run_id,library_name,abyss_version) values(null,'" << iV[2] << "','" << VERSION << "');";

			if (db.query (iStream.str()) && db.verbose_val > 1)
				std::cerr << iStream.str() << std::endl;
			std::ofstream ofile("db.txt");
			ofile << "done\n"; // so that iV.size != 3
		}
		std::stringstream uStream;
		uStream << "update Run_pe set stage=stage+1;";
		if (db.query (uStream.str()) && db.verbose_val >2)
			std::cerr << uStream.str() << std::endl;
		isrun = true; // as far as ifile exists
	}
	return isrun;
}

static inline std::string getPath (const std::string& program)
{
	std::stringstream echoStream, pathStream;
	echoStream << "echo `which " << program << "` > temp";
	std::string echoStr = echoStream.str();
	system(echoStr.c_str());
	std::ifstream ifile("temp");
	if (ifile) {
		std::string line;
		while (getline (ifile, line)) pathStream << line;
	}
	system("rm -f temp");
	return pathStream.str();
}

template <typename D>
static inline void assemblyStatsToDb (D& db, const std::string& program, const std::string& cmd)
{
	std::string temp = db.getProperTableName(program);
	std::stringstream sqlStream, mapKeys, mapValues;

	sqlStream << "\
create table if not exists " << temp << " (\
run_id integer default null, run_stage integer default null, exec_path text, command_line text, time_finish_" << temp << " not null default (datetime(current_timestamp,'localtime')), ";

	mapKeys << "run_id, run_stage, exec_path, command_line, ";
	mapValues << ((isRun(db)) ? "(select max(run_id) from Run_pe), (select stage from Run_pe where run_id = (select max(run_id) from Run_pe)), " : "null, null, ") << "'" << getPath(program) << "', '" << cmd << "', ";
	while (!dbMap.empty())
	{
		mapKeys << dbMap.begin()->first << (dbMap.size() == 1 ? "" : ", ");
		mapValues << dbMap.begin()->second << (dbMap.size() == 1 ? "" : ", ");
		sqlStream << dbMap.begin()->first << (dbMap.size() == 1 ? " int, foreign key(run_id) references Run_pe(run_id));" : " int, ");
		dbMap.erase(dbMap.begin());
	}

	// create tables
	if (db.query (sqlStream.str()) && db.verbose_val > 1)
		std::cerr << sqlStream.str() << std::endl;

	sqlStream.str("");
	sqlStream << "insert into " << temp << " (" << mapKeys.str() << ") values (" << mapValues.str() << ");";

	// insert values into tables
	if (db.query (sqlStream.str()) && db.verbose_val > 1)
		std::cerr << sqlStream.str() << std::endl;
}

template <typename Graph, typename D>
static inline void graphStatsToDb (const Graph& g, D& db, const std::string& program, const std::string& cmd)
{
	unsigned v = num_vertices(g) - num_vertices_removed(g);
	dbMap["V"] = v;

	unsigned e = num_edges(g);
	dbMap["E"] = e;

	assemblyStatsToDb(db, program, cmd);
}

static inline void addToDbMap (const std::string& key, const float& value)
{
	dbMap[key] = value;
}

template <typename Map>
static inline void addToDbMap (const Map& m)
{
	dbMap.insert (m.begin(), m.end());
}

#endif
