#ifndef WRITER_H
#define WRITER_H

#include <fstream>

class Writer
{
	public:
		Writer(const char* contigFilename);
		~Writer();
		
		// write out a contig
		void writeContig(const char* seq);
		
	private:
	
		int m_numContigs;
		
		// file handles
		std::fstream m_fContigs;
};

#endif
