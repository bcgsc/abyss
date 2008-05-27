#include <fstream>
#include "PairRecord.h"

PairRecord::PairRecord()
{

}

void PairRecord::addPairs(const PackedSeq& seq1, const PackedSeq& seq2)
{
	m_pairLookup[seq1].push_back(seq2);
	m_pairLookup[seq2].push_back(seq1);	
}

bool PairRecord::checkForPairs(const PackedSeq& seq) const
{
	if(m_pairLookup.count(seq) > 0)
	{
		return true;
	}
	else
	{
		return false;
	}
}

const PSequenceVector PairRecord::getPairs(const PackedSeq& seq) const
{
	PSequenceVector ret;
	addPairs(seq, ret);
	addPairs(reverseComplement(seq), ret);
	return ret;
}

void PairRecord::addPairs(const PackedSeq& seq, PSequenceVector& vec) const
{
	std::map<PackedSeq, PSequenceVector >::const_iterator findIter = m_pairLookup.find(seq);
	if(findIter != m_pairLookup.end())
	{
		for(PSequenceVector::const_iterator iter = findIter->second.begin(); iter != findIter->second.end(); ++iter)
		{
			vec.push_back(*iter);
		}
	}
	return;
}

void PairRecord::serialize(std::string filename)
{
	std::ofstream filehandle(filename.c_str(),std::ios::out | std::ios::binary);
	size_t numWrote = 0;
	for(PairLookupTable::const_iterator iter = m_pairLookup.begin(); iter != m_pairLookup.end(); ++iter)
	{
		// output the number of sequences
		size_t numSeqs = iter->second.size() + 1;
		
		filehandle.write((char*)&numSeqs, sizeof(numSeqs));
		
		// create a buffer
		size_t bufferSize = sizeof(PackedSeq) * numSeqs;
		char* buffer = new char[bufferSize];
		
		// write the key to the buffer
		size_t pos = 0;
		iter->first.serialize(buffer + pos);
		
		pos += sizeof(PackedSeq);
		// output its pairs
		for(PSequenceVector::const_iterator pairIter = iter->second.begin(); pairIter != iter->second.end(); ++pairIter)
		{
			pairIter->serialize(buffer + pos);
			pos += sizeof(PackedSeq);
		}
		
		assert(pos == bufferSize);
		
		// write the data
		filehandle.write(buffer, bufferSize);
		
		numWrote += numSeqs;
		
		// free the buffer
		delete [] buffer;
	}
	
	printf("wrote %zu sequences\n", numWrote);
	printf("num keys: %zu\n", m_pairLookup.size());	
	filehandle.close();
}

void PairRecord::unserialize(std::string filename)
{
	printf("Reading pairs from binary file %s\n", filename.c_str());
	std::ifstream filehandle(filename.c_str(), std::ios::binary);
	
	size_t numRead = 0;
	
	while(!filehandle.eof() && filehandle.peek() != EOF)
	{
		// read in the number of sequences in this record
		size_t numSeqs;
		filehandle.read((char*)&numSeqs, sizeof(numSeqs));
		
		// did we go past the end of the file
		if(filehandle.eof())
		{
			break;
		}
		
		// allocate space
		size_t bufferSize = sizeof(PackedSeq) * numSeqs;
		char* buffer = new char[bufferSize];
		
		// read
		filehandle.read(buffer, bufferSize);
		
		size_t pos = 0;
		// read the first (key) seq
		PackedSeq key;
		key.unserialize(buffer + pos);
		
		pos += sizeof(PackedSeq);
		
		PSequenceVector seqVec;
		seqVec.reserve(numSeqs - 1);
		
		for(size_t idx = 1; idx < numSeqs; ++idx)
		{
			// Could be faster
			PackedSeq curr;
			curr.unserialize(buffer + pos);
			seqVec.push_back(curr);
			
			pos += sizeof(PackedSeq);
		}
		
		m_pairLookup[key] = seqVec;
		
		numRead += numSeqs;
		
		delete [] buffer;		
	}
	
	printf("read %zu seqs\n", numRead);
	printf("num keys: %zu\n", m_pairLookup.size());
}
