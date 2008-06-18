#include "SetOperations.h"
#include "AlignmentCache.h"
#include <fstream>

AlignmentCache::AlignmentCache()
{
	
}

void AlignmentCache::addAlignment(const PackedSeq& seq, const ContigID& id, int position)
{
	//std::cout << "Adding " << seq.decode() << " " << id << " " << position << std::endl;
	AlignData alignment;
	alignment.contigID = id;
	alignment.position = position;
	alignment.isRC = false;
	
	m_alignmentCache[seq].insert(alignment);
	
	// insert the reverse complement as well
	alignment.isRC = true;
	m_alignmentCache[reverseComplement(seq)].insert(alignment);
}

void AlignmentCache::removeAlignment(const PackedSeq& seq, const ContigID& id, int position)
{
	// Remove the sequence and its reverse complement
	removeAlignmentInternal(seq, id, position);
	removeAlignmentInternal(reverseComplement(seq), id, position);
}

void AlignmentCache::removeAlignmentInternal(const PackedSeq& seq, const ContigID& id, int position)
{
	//std::cout << "Removing " << seq.decode() << " " << id << " " << position << std::endl;
	
	// Make sure this sequence has alignments
	AlignDB::iterator alignSetIter = m_alignmentCache.find(seq);
	assert(alignSetIter != m_alignmentCache.end());
	
	// Set up the align data structure
	AlignData lookupAlign;
	lookupAlign.contigID = id;
	lookupAlign.position = position;
	lookupAlign.isRC = false;
	
	// check if exists
	AlignSet::iterator matchIter = alignSetIter->second.find(lookupAlign);
	
	if(matchIter != alignSetIter->second.end())
	{
		alignSetIter->second.erase(matchIter);
	}
	else
	{
		std::cout << "Alignment not found! " << lookupAlign << " Alignments: " << std::endl;
		assert(false);
	}
}

void AlignmentCache::getAlignments(const PackedSeq& seq, AlignSet& outset) const
{
	// Make sure this sequence has alignments
	AlignDB::const_iterator alignDBIter = m_alignmentCache.find(seq);
	
	if(alignDBIter == m_alignmentCache.end())
	{
		printf("Get alignment failed for %s\n", seq.decode().c_str());
	}
	assert(alignDBIter != m_alignmentCache.end());
	
	outset = alignDBIter->second;
}

void AlignmentCache::translatePosition(const PackedSeq& seq, const ContigID& id, int oldPosition, int newPosition)
{
	//printf("Updating %s (%d->%d)\n", seq.decode().c_str(), oldPosition, newPosition);
	
	// Remove the old alignment
	removeAlignment(seq, id, oldPosition);
	
	// add the new alignment
	addAlignment(seq, id, newPosition);
}

bool AlignmentCache::compare(const AlignmentCache& otherDB)
{
	for(AlignDB::const_iterator keyIter = m_alignmentCache.begin(); keyIter != m_alignmentCache.end(); ++keyIter)
	{
		// make sure the key is in the other DB
		AlignDB::const_iterator dataIter = otherDB.m_alignmentCache.find(keyIter->first);
		if(dataIter == otherDB.m_alignmentCache.end())
		{
			printf("key not found! %s", keyIter->first.decode().c_str());
			return false;
		}
	}
	return true;
}

void AlignmentCache::printAlignmentsForSeq(const PackedSeq& seq) const
{
	AlignSet aSet;
	getAlignments(seq, aSet);
	
	for(AlignSet::iterator iter = aSet.begin(); iter != aSet.end(); ++iter)
	{
		std::cout << "Align: " << *iter << std::endl;
	}
}

void AlignmentCache::concatSets(ContigIDSet& seqSet1, const ContigIDSet& seqSet2) const
{
	seqSet1.insert(seqSet2.begin(), seqSet2.end());
}

void AlignmentCache::getSet(const PackedSeq& seq, ContigIDSet& outset) const
{
	(void)seq;
	(void)outset;
	/*
	AlignDB::const_iterator iter = m_alignmentCache.find(seq);
	if(iter != m_alignmentCache.end())
	{
		concatSets(outset, iter->second);
	}
	return;
	*/
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
