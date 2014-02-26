// DB_config.cc
#include <iostream>
#include <fstream>
#include <string>
using namespace std;

int main () {
	string species, strain, library;
	ofstream outputFile("db.txt");

	cout << "Enter a species name: ";
	cin >> species;

	cout << "Enter a strain name: ";
	cin >> strain;

	cout << "Enter a library name: ";
	cin >> library;

	outputFile << species << endl;
	outputFile << strain << endl;
	outputFile << library << endl;

	return 0;
}
