#include "Writer.h"

Writer::Writer(const char* contigFilename) : m_numContigs(0)
{
	m_fContigs.open(contigFilename, std::ios::out);
	assert(m_fContigs.is_open());
}

Writer::~Writer()
{
	m_fContigs.close();	
}

void Writer::writeContig(const char* seq)
{
	char header[64];
	sprintf(header, ">contig_%d %d\n", m_numContigs, strlen(seq));
	
	m_fContigs << header;
	m_fContigs << seq << "\n";
	m_fContigs.flush();
	m_numContigs++;
}
