#ifndef CONFIG_H
#define CONFIG_H

#include <string>

class Config
{
	public:
		Config();
		void readConfig(const std::string filename);
		
		std::string getRootDataDir() const;
		int getSequenceLength() const;
		int getPartitionStep() const;
	
	private:		
		void parseNameValue(const char* buffer, std::string& name, std::string& value);
		
		std::string m_rootDataDirectory;
		int m_partitionStep;
		int m_sequenceLength;
};

#endif
