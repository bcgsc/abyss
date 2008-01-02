#ifndef CONFIG_H
#define CONFIG_H

#include <string>

class Config
{
	public:
		Config();
		void readConfig(const std::string filename);
		
		std::string getRootDataDir() const;
		std::string getSequenceFilename() const;
		int getSequenceLength() const;
		int getUnitSize() const;
	
	private:		
		void parseNameValue(const char* buffer, std::string& name, std::string& value);
		
		std::string m_rootDataDirectory;
		std::string m_sequenceFilename;
		int m_unitSize;
		int m_sequenceLength;
};

#endif
