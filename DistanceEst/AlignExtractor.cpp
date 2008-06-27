#include "AlignExtractor.h"

AlignExtractor::AlignExtractor(std::string file)
{
	m_fileHandle.open(file.c_str());

	// Prime the read by reading in the first contig
	m_currPair = readRecord();
}

AlignExtractor::~AlignExtractor()
{
	m_fileHandle.close();
}

bool AlignExtractor::extractContigAlignments(AlignPairVec& outPairs)
{
	AlignPair previousPair = m_currPair;

	bool isEOF = false;
	while(previousPair.refRec.contig == m_currPair.refRec.contig && !isEOF)
	{
		outPairs.push_back(m_currPair);
		m_currPair = readRecord();
		
		isEOF = (m_fileHandle.eof() || m_fileHandle.peek() == EOF);
	}
	return isEOF;
}

// Read a single record in
AlignPair AlignExtractor::readRecord()
{

	std::string readID1;
	std::string readID2;
	/*
	std::string line;
	getline(m_fileHandle, line);
	std::cout << "LINE: " << line << std::endl;
	*/
	Alignment ali1;
	Alignment ali2;

	m_fileHandle >> readID1 >> ali1 >> readID2 >> ali2;

	AlignPair ap = {ali1, ali2};
	return ap;
}
