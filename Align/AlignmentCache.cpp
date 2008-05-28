#include "AlignmentCache.h"
#include <fstream>

AlignmentCache::AlignmentCache()
{
	
}

void AlignmentCache::addAlignment(const PackedSeq& seq, const ContigID& id)
{
	m_alignmentCache[seq].insert(id);
}

/*
void AlignmentCache::serialize(std::string filename)
{
	std::ofstream filehandle(filename.c_str(),std::ios::out | std::ios::binary);
	size_t numWrote = 0;
	for(AlignDB::const_iterator iter = m_alignmentCache.begin(); iter != m_alignmentCache.end(); ++iter)
	{
		// write the key to the file
		const size_t keySize = sizeof(PackedSeq);
		char keyBuffer[keySize];
		iter->first.serialize(keyBuffer);
		filehandle.write(keyBuffer, keySize);
		
		const size_t numAligns = iter->second.size();
		
		// write the number of alignments to the file
		filehandle.write((char*)&numAligns, sizeof(numAligns));
		
		// create a buffer
		const size_t alignSize = sizeof(ContigID);
		size_t bufferSize = alignSize * numAligns;
		char* buffer = new char[bufferSize];
		size_t pos = 0;
		for(ContigIDColl::const_iterator alignIter = iter->second.begin(); alignIter != iter->second.end(); ++alignIter)
		{
			// write the alignments to the buffer
			memcpy(buffer + pos, &*alignIter, alignSize);
			pos += alignSize;
		}
				
		assert(pos == bufferSize);
		
		// write the data
		filehandle.write(buffer, bufferSize);
		
		numWrote += numAligns;
		
		// free the buffer
		delete [] buffer;
	}
	
	printf("wrote %zu alignments\n", numWrote);
	printf("num keys: %zu\n", m_alignmentCache.size());	
	filehandle.close();
}

void AlignmentCache::unserialize(std::string filename)
{
	printf("Reading alignments from binary file %s\n", filename.c_str());
	std::ifstream filehandle(filename.c_str(), std::ios::binary);
	
	size_t numRead = 0;
	
	while(!filehandle.eof() && filehandle.peek() != EOF)
	{
		// read in the key
		PackedSeq key;
		const size_t keySize = sizeof(key);
		char keyBuffer[keySize];
		
		filehandle.read(keyBuffer, keySize);
		key.unserialize(keyBuffer);
		
		// did we go past the end of the file
		if(filehandle.eof())
		{
			break;
		}
				
		// read in the number of sequences in this record
		size_t numAlignments;
		filehandle.read((char*)&numAlignments, sizeof(numAlignments));
		

		// allocate space
		const size_t alignSize = sizeof(Alignment);
		size_t bufferSize = alignSize * numAlignments;
		char* buffer = new char[bufferSize];
		
		// read
		filehandle.read(buffer, bufferSize);
		
		size_t pos = 0;
		
		// copy the alignments out
		std::vector aligns;
		aligns.reserve(numAlignments);

		for(size_t idx = 0; idx < numAlignments; ++idx)
		{
			Alignment curr;
			memcpy(&curr, buffer + pos, alignSize);
			aligns.push_back(curr);
			
			pos += alignSize;
		}
		
		m_alignmentCache[key].insert;
		
		numRead += numAlignments;
		
		delete [] buffer;		
	}
	
	printf("read %zu alignments\n", numRead);
	printf("num keys: %zu\n", m_alignmentCache.size());
}
*/
