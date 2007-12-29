#include <fstream>
#include "Config.h"

#define PARTITION_DIR_NAME "partition_dir"
#define PARTITION_STEP "partition_step"
#define SEQUENCE_LENGTH "sequence_length"

Config::Config() : m_partitionStep(-1), m_sequenceLength(-1)
{
	
}

void Config::readConfig(const std::string filename)
{
	std::ifstream inFile(filename.c_str());
	assert(inFile.is_open());
	
	const int MAX_LINE_SIZE = 512;
	char buffer[MAX_LINE_SIZE];
	
	while(!inFile.eof() && inFile.peek() != EOF)
	{
		// Parse the config value	
		inFile.getline(buffer, MAX_LINE_SIZE);
		
		// each line should be a name/value pairs seperated by "="
		std::string name;
		std::string value;
		parseNameValue(buffer, name, value);
		
		if(name == PARTITION_DIR_NAME)
		{
			m_rootDataDirectory = value;
		}
		else if(name == PARTITION_STEP)
		{
			m_partitionStep = atoi(value.c_str());
		}
		else if(name == SEQUENCE_LENGTH)
		{
			m_sequenceLength = atoi(value.c_str());
		}
		else
		{
			printf("could not parse config name: %s value: %s\n", name.c_str(), value.c_str());
			assert(false);
		}	
	}
}

void Config::parseNameValue(const char* buffer, std::string& name, std::string& value)
{
	std::string inStr(buffer);
	int pos = inStr.find('=');
	
	name = inStr.substr(0, pos);
	value = inStr.substr(pos+1, inStr.size() - pos);
}

std::string Config::getRootDataDir() const
{
	return m_rootDataDirectory;
}

int Config::getSequenceLength() const
{
	return m_sequenceLength;
}

int Config::getPartitionStep() const
{
	return m_partitionStep;
}
