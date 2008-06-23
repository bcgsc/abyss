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
	while(previousPair.refRec.contigID == m_currPair.refRec.contigID && !isEOF)
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
	ContigID id1;
	ContigID id2;
	std::string readID1;
	std::string readID2;
	int pos1;
	int pos2;
	bool rc1;
	bool rc2;
	std::string discard;

	// ref read
	m_fileHandle >> readID1 >> discard >> id1 >> discard >> pos1 >> discard >> rc1;

	// pair read
	m_fileHandle >> readID2 >> discard >> id2 >> discard >> pos2 >> discard >> rc2;

	// build aligns
	AlignData ref = {id1, pos1, rc1};
	AlignData pair = {id2, pos2, rc2};

	AlignPair ap = {ref, pair};
	return ap;
}
